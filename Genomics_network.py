# -*- coding: utf-8 -*-
"""
Genomics — Network Propagation / Module Detection (STRING).

Unlike ORA/Burden (category-based: test genes against predefined gene sets),
this analysis is network-based: map the user's genes onto the STRING
protein-protein interaction network for the chosen organism, propagate each
seed's signal across the graph (random walk with restart), assess which
proteins accumulate more signal than a degree-preserving permutation null
would predict ("hits"), then optionally cluster the hits into densely
interconnected modules (Markov Clustering / MCL) and label each module by
functional enrichment (reusing the ORA Fisher/BH engine in Genomics_module.py).

Data source: STRING (string-db.org) bulk per-organism download files, CC BY
4.0 (free academic + commercial use, attribution only; redistribution and
hosting are permitted - verified against string-db.org/cgi/access). The live
STRING REST API is deliberately NOT used here: STRING's own documentation
asks that the API be reserved for "limited or occasional" access and directs
bulk/whole-organism use to the Downloads page instead - exactly what this
module does, mirroring the existing Reactome/GWAS bulk-download-and-cache
pattern in Genomics_module.py.

Four STRING files are used, all under https://stringdb-downloads.org/download/:
  - species.v12.0.txt                        - organism coverage list (drives
                                                the organism selector for this
                                                module specifically, since
                                                STRING covers a SUBSET of the
                                                taxa UniProt/GO/Reactome do).
  - {taxon}.protein.links.v12.0.txt.gz        - "protein1 protein2 combined_score"
                                                (space-separated; combined_score
                                                0-1000) -> the network edges.
  - {taxon}.protein.info.v12.0.txt.gz         - STRING protein id -> gene symbol.
  - {taxon}.protein.aliases.v12.0.txt.gz      - STRING protein id -> external
                                                identifiers; filtered to
                                                source == "UniProt_AC" gives
                                                the UniProt bridge used
                                                everywhere else in Omiscera.

Formats verified by downloading and inspecting real STRING v12.0 files during
development (not guessed from documentation alone).
"""

import os
import time
from collections import OrderedDict

import numpy as np
import pandas as pd

try:
    import scipy.sparse as sp
except ImportError:
    sp = None

import Genomics_module as genomics_ora

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

STRING_VERSION = "12.0"
STRING_DOWNLOAD_BASE = "https://stringdb-downloads.org/download"
STRING_SPECIES_URL = f"{STRING_DOWNLOAD_BASE}/species.v{STRING_VERSION}.txt"
# STRING releases roughly yearly, far less often than GO/Reactome/GWAS - a
# longer cache TTL than genomics_ora's 30-day defaults is appropriate.
STRING_REFRESH_DAYS = 90

STRING_CACHE_DIR = os.path.join(genomics_ora.CACHE_DIR, "string")
STRING_SPECIES_PATH = os.path.join(STRING_CACHE_DIR, f"species.v{STRING_VERSION}.txt")

DEFAULT_RESTART_PROB = 0.5       # random-walk-with-restart alpha (fixed, reported)
DEFAULT_MCL_INFLATION = 2.0      # MCL inflation parameter (fixed, reported)
MCL_EXPANSION = 2                # MCL expansion power (standard value; fixed)
DEFAULT_MIN_MODULE_SIZE = 3
DEFAULT_CONFIDENCE = 700         # STRING combined_score threshold, 0-1000 scale (== confidence 0.7)
CONFIDENCE_CHOICES = (150, 400, 700, 900)
DEFAULT_N_PERMUTATIONS = 1000
N_PERMUTATIONS_CHOICES = (500, 1000, 2000, 5000, 10000)
DEFAULT_FDR_ALPHA = 0.05
# Caps per-module enrichment term size so a module's label isn't dominated by
# a huge, uninformative catch-all GO term (e.g. "protein binding", ~10k
# genes) - the same practice standard enrichment tools (g:Profiler, Enrichr,
# DAVID) apply by default. Raised from 500 to 2000 once GO ancestor
# propagation (genomics_ora.build_go_universe) was added: propagation pushes
# genuinely useful, recognizable mid-level labels ("DNA repair", "cell
# cycle", "immune response") above the old 500-gene bound, so a 500 cap would
# now silently exclude exactly the labels a module-detection user wants and
# leave obscure deep leaf terms instead. 2000 still excludes the true
# catch-all/near-root terms. Fixed, but reported for reproducibility.
DEFAULT_MODULE_MAX_TERM_SIZE = 2000

RWR_MAX_ITER = 200
RWR_TOL = 1e-8
MCL_MAX_ITER = 100
MCL_PRUNE_THRESHOLD = 1e-5
MCL_CONVERGE_TOL = 1e-6
DEGREE_BIN_COUNT = 10             # for degree-preserving permutation seed resampling

NETPROP_HIT_COLUMNS = [
    "STRING ID", "Gene Symbol", "UniProt ID", "Propagated score",
    "Z-score", "p-value", "FDR", "Seed", "Seed score",
]
NETPROP_MODULE_COLUMNS = [
    "Module ID", "Size", "Genes", "Top enriched term", "Top enriched FDR", "Enrichment database",
]

# In-process cache of built (taxon, confidence) -> parsed network dict, so a
# second run in the same server process doesn't re-parse the ~14M-row links
# file (mirrors the retriever's in-memory BLAST_CACHE pattern; the raw
# downloaded files are ALSO cached on disk via _ensure_string_file below, so
# even a fresh process only re-downloads, never re-fetches from STRING).
# Bounded: a human graph is ~50-100 MB and this used to be an unbounded dict,
# so a long-lived server that saw several organisms (or several confidence
# cutoffs for one organism) grew without limit. Two entries is enough for the
# realistic case - re-running the same organism, or comparing two confidence
# settings - and eviction only costs a re-parse of already-downloaded files.
_GRAPH_CACHE: "OrderedDict[tuple, dict]" = OrderedDict()
_GRAPH_CACHE_MAX_ENTRIES = 2


def _graph_cache_put(key: tuple, value: dict) -> None:
    """Insert with LRU eviction. Kept tiny and separate so the eviction policy
    is visible and testable rather than buried in build_string_network."""
    _GRAPH_CACHE[key] = value
    _GRAPH_CACHE.move_to_end(key)
    while len(_GRAPH_CACHE) > _GRAPH_CACHE_MAX_ENTRIES:
        _GRAPH_CACHE.popitem(last=False)


# ---------------------------------------------------------------------------
# STRING data layer
# ---------------------------------------------------------------------------

def _ensure_string_file(url: str, local_path: str, refresh_days: float = STRING_REFRESH_DAYS, verbose: bool = False) -> str:
    """Download (if needed) and cache one STRING bulk file. Mirrors
    Genomics_module._ensure_reactome_file's download/cache pattern."""
    if os.path.exists(local_path):
        age_days = (time.time() - os.path.getmtime(local_path)) / 86400
        if age_days < refresh_days:
            if verbose:
                print(f"[STRING] Using cached {os.path.basename(local_path)} ({age_days:.1f} days old)")
            return local_path
        if verbose:
            print(f"[STRING] Cache expired for {os.path.basename(local_path)} ({age_days:.1f} days), re-downloading...")

    os.makedirs(os.path.dirname(local_path), exist_ok=True)
    assert genomics_ora.HTTP_SESSION is not None, "requests not installed"
    if verbose:
        print(f"[STRING] Downloading {url} ...")
    resp = genomics_ora.HTTP_SESSION.get(url, headers=genomics_ora.OMISCERA_HEADERS, timeout=180, stream=True)
    resp.raise_for_status()
    downloaded = 0
    with open(local_path, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=1024 * 256):
            fh.write(chunk)
            downloaded += len(chunk)
    if verbose:
        print(f"[STRING] Downloaded {downloaded / 1e6:.1f} MB -> {local_path}")
    return local_path


def load_string_species(verbose: bool = False) -> pd.DataFrame:
    """Load STRING's covered-organism list. Drives the organism selector for
    THIS module specifically: STRING covers a subset of the taxa UniProt/GO/
    Reactome cover, so sourcing the menu from STRING (rather than UniProt's
    taxonomy, as ORA/Burden do) guarantees a network exists for whatever the
    user picks. Columns (STRING v12.0): taxon_id, STRING_type,
    STRING_name_compact, official_name_NCBI, domain.
    """
    path = _ensure_string_file(STRING_SPECIES_URL, STRING_SPECIES_PATH, verbose=verbose)
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = [c.lstrip("#").strip() for c in df.columns]
    return df


def search_string_species(query: str, limit: int = 15, verbose: bool = False) -> list:
    """Case-insensitive substring search over STRING's species list, by taxon
    id or name. Returns [{"taxon_id": str, "name": str}, ...]."""
    q = str(query).strip().lower()
    if not q:
        return []
    df = load_string_species(verbose=verbose)
    name_col = "official_name_NCBI" if "official_name_NCBI" in df.columns else "STRING_name_compact"
    mask = (
        df["taxon_id"].astype(str).str.lower().str.contains(q, na=False, regex=False)
        | df[name_col].astype(str).str.lower().str.contains(q, na=False, regex=False)
    )
    hits = df[mask].head(limit)
    return [
        {"taxon_id": str(row["taxon_id"]), "name": str(row.get(name_col, ""))}
        for _, row in hits.iterrows()
    ]


def is_taxon_in_string(taxon_id: str, verbose: bool = False) -> bool:
    """True if STRING has a network for this NCBI taxon id."""
    df = load_string_species(verbose=verbose)
    return str(taxon_id) in set(df["taxon_id"].astype(str))


def _string_file_url(file_type: str, taxon_id: str) -> str:
    return f"{STRING_DOWNLOAD_BASE}/{file_type}.v{STRING_VERSION}/{taxon_id}.{file_type}.v{STRING_VERSION}.txt.gz"


def _string_file_paths(taxon_id: str) -> dict:
    return {
        "links": os.path.join(STRING_CACHE_DIR, f"{taxon_id}.protein.links.v{STRING_VERSION}.txt.gz"),
        "info": os.path.join(STRING_CACHE_DIR, f"{taxon_id}.protein.info.v{STRING_VERSION}.txt.gz"),
        "aliases": os.path.join(STRING_CACHE_DIR, f"{taxon_id}.protein.aliases.v{STRING_VERSION}.txt.gz"),
    }


def build_string_network(
    taxon_id: str,
    confidence: int = DEFAULT_CONFIDENCE,
    verbose: bool = False,
    force_refresh: bool = False,
) -> dict:
    """Download (if needed), cache, and parse the STRING network for one
    organism at a given confidence cutoff (STRING combined_score, 0-1000).

    Returns:
      "adjacency":         scipy.sparse.csr_matrix, column-normalized (each
                            column sums to 1) - the transition matrix used by
                            the random walk.
      "raw_adjacency":      scipy.sparse.csr_matrix, symmetric, weights =
                            combined_score / 1000 - used for degree binning
                            and as the module-detection input (unnormalized
                            confidence weights).
      "protein_index":      {string_protein_id: matrix row/col index}.
      "index_protein":      [string_protein_id, ...] (row i -> id).
      "protein_gene":       {string_protein_id: preferred_name}.
      "uniprot_to_protein": {UniProt accession: string_protein_id} (from the
                            aliases file's UniProt_AC rows; a STRING protein
                            can have multiple UniProt accessions, e.g.
                            isoforms - last one wins, an accepted
                            approximation mirroring other identifier bridges
                            in this codebase).
      "n_proteins":         int.
    """
    cache_key = (str(taxon_id), int(confidence))
    if not force_refresh and cache_key in _GRAPH_CACHE:
        _GRAPH_CACHE.move_to_end(cache_key)
        return _GRAPH_CACHE[cache_key]
    if sp is None:
        raise RuntimeError("scipy is required for network propagation")

    paths = _string_file_paths(taxon_id)
    links_path = _ensure_string_file(_string_file_url("protein.links", taxon_id), paths["links"], verbose=verbose)
    info_path = _ensure_string_file(_string_file_url("protein.info", taxon_id), paths["info"], verbose=verbose)
    aliases_path = _ensure_string_file(_string_file_url("protein.aliases", taxon_id), paths["aliases"], verbose=verbose)

    if verbose:
        print(f"[STRING] Parsing links for taxon {taxon_id} at confidence >= {confidence}...")
    result = _parse_string_network_files(links_path, info_path, aliases_path, confidence, taxon_id=taxon_id)
    _graph_cache_put(cache_key, result)
    if verbose:
        print(f"[STRING] Built network: {result['n_proteins']} proteins (>= {confidence})")
    return result


def _parse_string_network_files(links_path: str, info_path: str, aliases_path: str, confidence: int, taxon_id: str = "") -> dict:
    """Parse already-downloaded STRING bulk files into the graph dict
    build_string_network returns. Extracted from build_string_network
    (identical logic, just given a name) so the parsing/symmetrization/
    normalization/alias-resolution logic is independently unit-testable on
    small hand-built fixture files, without a real multi-hundred-MB download.
    See build_string_network's docstring for the return shape.
    """
    links_df = pd.read_csv(links_path, sep=" ", dtype={"protein1": str, "protein2": str, "combined_score": "int32"})
    links_df = links_df[links_df["combined_score"] >= confidence]
    if links_df.empty:
        raise ValueError(f"No STRING interactions for taxon {taxon_id} at confidence >= {confidence}.")

    proteins = pd.unique(pd.concat([links_df["protein1"], links_df["protein2"]], ignore_index=True))
    protein_index = {p: i for i, p in enumerate(proteins)}
    n = len(proteins)

    rows = links_df["protein1"].map(protein_index).to_numpy()
    cols = links_df["protein2"].map(protein_index).to_numpy()
    weights = links_df["combined_score"].to_numpy(dtype=np.float64) / 1000.0

    raw_adj = sp.csr_matrix((weights, (rows, cols)), shape=(n, n))
    raw_adj = raw_adj.maximum(raw_adj.T)  # symmetrize defensively (STRING already lists both directions)

    # Column-normalize for the random walk: p_{t+1} = W @ p_t must redistribute
    # probability mass without gaining/losing total mass, so each column of W
    # sums to 1. A zero-degree column (isolated node post-filtering) is left
    # as an all-zero column - it neither sends nor receives signal, correct.
    col_sums = np.asarray(raw_adj.sum(axis=0)).flatten()
    col_sums_safe = np.where(col_sums == 0, 1.0, col_sums)
    inv_diag = sp.diags(1.0 / col_sums_safe)
    # sp.csr_matrix(...) (not .tocsr() on the matmul result) - scipy's sparse
    # __matmul__ stubs return an imprecise union type Pylance can't narrow to
    # something with .tocsr(); wrapping in the constructor is equivalent at
    # runtime (an already-CSR matrix passes through) and resolves cleanly.
    norm_adj = sp.csr_matrix(raw_adj @ inv_diag)

    info_df = pd.read_csv(info_path, sep="\t", dtype=str)
    info_df.columns = [c.lstrip("#").strip() for c in info_df.columns]
    protein_gene = dict(zip(info_df["string_protein_id"], info_df["preferred_name"]))

    aliases_df = pd.read_csv(aliases_path, sep="\t", dtype=str)
    aliases_df.columns = [c.lstrip("#").strip() for c in aliases_df.columns]
    uniprot_rows = aliases_df[aliases_df["source"] == "UniProt_AC"]
    uniprot_to_protein = dict(zip(uniprot_rows["alias"], uniprot_rows["string_protein_id"]))

    return {
        "adjacency": norm_adj,
        "raw_adjacency": raw_adj,
        "protein_index": protein_index,
        "index_protein": list(proteins),
        "protein_gene": protein_gene,
        "uniprot_to_protein": uniprot_to_protein,
        "n_proteins": n,
    }


# ---------------------------------------------------------------------------
# Random walk with restart (propagation)
# ---------------------------------------------------------------------------

def propagate_rwr(
    adjacency,
    seed_vector: np.ndarray,
    restart_prob: float = DEFAULT_RESTART_PROB,
    max_iter: int = RWR_MAX_ITER,
    tol: float = RWR_TOL,
) -> np.ndarray:
    """Random walk with restart: p_{t+1} = (1-alpha)*W@p_t + alpha*s, where W
    is the column-normalized adjacency, s is the seed vector normalized to
    sum to 1, and alpha = restart_prob (higher alpha -> signal stays closer to
    the seeds; lower alpha -> spreads further through the network). Iterates
    to L1 convergence or max_iter. Returns the converged probability vector
    (sums to 1); a protein's value is its accumulated propagated signal.
    """
    s = np.asarray(seed_vector, dtype=np.float64)
    total = s.sum()
    if total <= 0:
        raise ValueError("seed_vector must have at least one positive entry")
    s = s / total
    p = s.copy()
    alpha = restart_prob
    for _ in range(max_iter):
        p_next = (1 - alpha) * (adjacency @ p) + alpha * s
        if np.abs(p_next - p).sum() < tol:
            p = p_next
            break
        p = p_next
    return p


def _degree_bins(raw_adjacency, n_bins: int = DEGREE_BIN_COUNT):
    """Bin protein indices into n_bins degree-quantile bins (by summed raw
    edge weight), for degree-preserving random seed resampling in the
    permutation null below. Returns (bin_of_index_array, {bin: [indices]})."""
    degrees = np.asarray(raw_adjacency.sum(axis=0)).flatten()
    n = len(degrees)
    order = np.argsort(degrees)
    bin_of = np.zeros(n, dtype=int)
    bin_size = max(1, n // n_bins)
    for b in range(n_bins):
        start = b * bin_size
        end = (b + 1) * bin_size if b < n_bins - 1 else n
        bin_of[order[start:end]] = b
    bins_map: dict = {}
    for idx, b in enumerate(bin_of):
        bins_map.setdefault(int(b), []).append(idx)
    return bin_of, bins_map


def _welford_update(mean, m2, count, x):
    """One step of Welford's online mean/variance update, extracted from the
    permutation-null loop below so it's independently unit-testable. count
    is the UPDATED (1-indexed) sample count after including x; the caller
    reads variance back out as m2 / (count - 1)."""
    delta = x - mean
    mean = mean + delta / count
    m2 = m2 + delta * (x - mean)
    return mean, m2


def _map_seeds_to_indices(id_map: dict, protein_index: dict, seed_scores: dict) -> tuple:
    """Map external identifiers onto graph row indices.

    id_map: {external id: string_protein_id}. Genomics passes
    graph["uniprot_to_protein"]; the Transcriptomics network module passes a
    gene-symbol (or Ensembl/Entrez) bridge instead. Extracting this is what
    makes the numerical core below identifier-agnostic - nothing downstream of
    here knows or cares which identifier space the caller used.

    Two behaviours are load-bearing and preserved exactly:
      - several external ids collapsing onto ONE protein have their scores
        SUMMED, not overwritten;
      - if any resulting weight is <= 0 the whole set is shifted positive by
        |min| + 1, because RWR redistributes probability mass and needs
        strictly positive seeds.

    That shift is monotonic but NOT scale-preserving: propagate_rwr normalizes
    the seed vector to sum 1, so only weight RATIOS survive, and shifting
    changes them. A caller passing signed values (e.g. a signed ranking metric)
    therefore silently gets a different analysis than it intended. Callers that
    guarantee positive weights never reach the shift; the Transcriptomics module
    guarantees exactly that, deliberately.

    Returns (seed_idx_score {index: positive weight}, n_not_mapped).
    """
    seed_idx_score: dict = {}
    n_not_mapped = 0
    for ext_id, score in seed_scores.items():
        pid = id_map.get(ext_id)
        if pid is None or pid not in protein_index:
            n_not_mapped += 1
            continue
        idx = protein_index[pid]
        seed_idx_score[idx] = seed_idx_score.get(idx, 0.0) + float(score)

    if seed_idx_score:
        min_score = min(seed_idx_score.values())
        if min_score <= 0:
            shift = abs(min_score) + 1.0
            seed_idx_score = {k: v + shift for k, v in seed_idx_score.items()}
    return seed_idx_score, n_not_mapped


def _propagation_core(
    graph: dict,
    seed_idx_score: dict,
    restart_prob: float = DEFAULT_RESTART_PROB,
    n_permutations: int = DEFAULT_N_PERMUTATIONS,
    seed_rng: int = 42,
    progress_callback=None,
) -> dict:
    """The numerical heart: observed RWR, degree-preserving permutation null,
    z-scores, normal-tail p-values, BH-FDR. Works purely in index space and
    contains NO identifier vocabulary, so both the Genomics (UniProt) and
    Transcriptomics (gene symbol) modules share one implementation rather than
    keeping two copies that could drift apart and give different answers to the
    same scientific question.

    Standardizing against a per-protein null (rather than counting how often a
    random seed set exceeds the observed) is deliberate: on a dense PPI network
    hubs soak up signal from ANY high-degree seeds, so a raw exceedance count is
    far too conservative. Z-scoring each protein against its own null removes
    that baseline. Welford's online mean/variance keeps memory O(n).

    Returns observed/zscores/pvalues/fdrs plus two honesty diagnostics:
    n_exceed (how many permutations actually met or beat the observed - the
    empirical view the z-scores approximate) and null_skew_at_max (skewness of
    the null at the top-scoring protein, since the p-values assume a normal
    null that a permutation distribution is typically right-skewed against).
    Neither is used as a p-value: at ~1000 permutations the empirical floor of
    1/(n+1), after BH across ~20k proteins, would call nothing significant.
    """
    def _progress(current, total, label):
        if progress_callback:
            try:
                progress_callback(current, total, label)
            except Exception:
                pass

    adjacency = graph["adjacency"]
    raw_adjacency = graph["raw_adjacency"]
    n = graph["n_proteins"]

    seed_vector = np.zeros(n, dtype=np.float64)
    for idx, score in seed_idx_score.items():
        seed_vector[idx] = score

    _progress(0, n_permutations + 1, "Propagating seed signal...")
    observed = propagate_rwr(adjacency, seed_vector, restart_prob=restart_prob)

    bin_of, bins_map = _degree_bins(raw_adjacency)
    seed_indices = list(seed_idx_score.keys())
    seed_scores_list = [seed_idx_score[i] for i in seed_indices]

    rng = np.random.default_rng(seed_rng)
    null_mean = np.zeros(n, dtype=np.float64)
    null_m2 = np.zeros(n, dtype=np.float64)
    n_exceed = np.zeros(n, dtype=np.int64)
    top_idx = int(np.argmax(observed))
    top_null_samples = []
    count = 0
    for perm_i in range(n_permutations):
        perm_vector = np.zeros(n, dtype=np.float64)
        used: set = set()
        for orig_idx, score in zip(seed_indices, seed_scores_list):
            candidates = bins_map.get(int(bin_of[orig_idx]), list(range(n)))
            cand = candidates[rng.integers(0, len(candidates))]
            for _try in range(10):
                if cand not in used:
                    break
                cand = candidates[rng.integers(0, len(candidates))]
            used.add(cand)
            perm_vector[cand] += score
        perm_p = propagate_rwr(adjacency, perm_vector, restart_prob=restart_prob)
        count += 1
        null_mean, null_m2 = _welford_update(null_mean, null_m2, count, perm_p)
        n_exceed += (perm_p >= observed)
        top_null_samples.append(float(perm_p[top_idx]))
        if (perm_i + 1) % max(1, n_permutations // 20) == 0:
            _progress(perm_i + 1, n_permutations + 1, f"Permutation {perm_i + 1}/{n_permutations}")

    null_std = np.sqrt(null_m2 / max(count - 1, 1))
    with np.errstate(divide="ignore", invalid="ignore"):
        zscores = (observed - null_mean) / null_std
    # Proteins with zero null variability (e.g. never reached by any random
    # seed set) are treated as non-significant rather than +inf - a
    # conservative default that avoids spurious hits from numerical edge cases.
    zscores = np.where(null_std > 0, zscores, 0.0)
    # One-sided upper-tail p-value from the standard normal (ndtr is the normal
    # CDF; it's a vectorized special function, not a BLAS matmul, so it is not
    # affected by the dense-BLAS issue noted in mcl_cluster).
    from scipy.special import ndtr
    pvalues = ndtr(-zscores)
    fdrs = genomics_ora.benjamini_hochberg(list(pvalues))

    skew = 0.0
    if len(top_null_samples) > 2:
        arr = np.asarray(top_null_samples, dtype=np.float64)
        sd = arr.std(ddof=1)
        if sd > 0:
            skew = float((((arr - arr.mean()) / sd) ** 3).mean())

    return {
        "observed": observed,
        "null_mean": null_mean,
        "null_std": null_std,
        "zscores": zscores,
        "pvalues": pvalues,
        "fdrs": fdrs,
        "n_exceed": n_exceed,
        "n_exceed_top": int(n_exceed[top_idx]),
        "null_skew_at_max": skew,
        "n_permutations": count,
        "seed_idx_score": seed_idx_score,
    }

def run_propagation_with_significance(
    graph: dict,
    seed_uniprot_scores: dict,
    restart_prob: float = DEFAULT_RESTART_PROB,
    n_permutations: int = DEFAULT_N_PERMUTATIONS,
    fdr_alpha: float = DEFAULT_FDR_ALPHA,
    seed_rng: int = 42,
    progress_callback=None,
    verbose: bool = False,
) -> dict:
    """Map seeds onto the STRING network, run RWR, and assess each protein's
    propagated score against a degree-preserving permutation null (shuffle
    WHICH proteins are seeds - matched by degree bin, so a hub doesn't look
    "significant" merely for being a hub - while reusing the same seed
    weights). This is the standard network-propagation significance approach
    (cf. HotNet2-style permutation testing).

    seed_uniprot_scores: {UniProt accession: score}. Scores must already be
        ORIENTED so higher = more signal (the caller negates upstream for a
        "lower = more important" metric, same convention as the Burden Test).
        Pass all 1.0 for binary/equal-weight seeds (no score column supplied).

    Returns a dict with "table" (NETPROP_HIT_COLUMNS, every scored protein -
    callers filter to FDR < fdr_alpha for "hits"), "n_seeds_mapped",
    "n_seeds_not_in_string", "n_hits", "n_recruited", and "seed_string_ids"
    (set of STRING ids that were seeds - feeds module detection).
    """
    def _progress(current, total, label):
        if progress_callback:
            try:
                progress_callback(current, total, label)
            except Exception:
                pass

    # adjacency / raw_adjacency are read inside _propagation_core now; this
    # wrapper only handles identifier vocabulary and table assembly.
    protein_index = graph["protein_index"]
    index_protein = graph["index_protein"]
    uniprot_to_protein = graph["uniprot_to_protein"]
    protein_gene = graph["protein_gene"]
    n = graph["n_proteins"]

    string_to_uniprot = {v: k for k, v in uniprot_to_protein.items()}

    seed_idx_score, n_not_in_string = _map_seeds_to_indices(
        uniprot_to_protein, protein_index, seed_uniprot_scores)
    n_seeds_mapped = len(seed_idx_score)
    if n_seeds_mapped == 0:
        raise ValueError("None of the input genes mapped to a STRING protein for this organism/confidence.")

    core = _propagation_core(
        graph, seed_idx_score, restart_prob=restart_prob, n_permutations=n_permutations,
        seed_rng=seed_rng, progress_callback=progress_callback)
    observed = core["observed"]
    zscores = core["zscores"]
    pvalues = core["pvalues"]
    fdrs = core["fdrs"]
    seed_idx_score = core["seed_idx_score"]

    rows = []
    for i in range(n):
        pid = index_protein[i]
        is_seed = i in seed_idx_score
        rows.append({
            "STRING ID": pid,
            "Gene Symbol": protein_gene.get(pid, ""),
            "UniProt ID": string_to_uniprot.get(pid, ""),
            "Propagated score": float(observed[i]),
            "Z-score": float(zscores[i]),
            "p-value": float(pvalues[i]),
            "FDR": float(fdrs[i]),
            "Seed": bool(is_seed),
            "Seed score": float(seed_idx_score[i]) if is_seed else None,
        })
    table = pd.DataFrame(rows, columns=NETPROP_HIT_COLUMNS)
    table = table.sort_values(by="FDR", ascending=True, kind="stable").reset_index(drop=True)

    hits_mask = table["FDR"] < fdr_alpha
    n_hits = int(hits_mask.sum())
    n_recruited = int((hits_mask & ~table["Seed"]).sum())
    seed_string_ids = {index_protein[i] for i in seed_idx_score}

    _progress(n_permutations + 1, n_permutations + 1, "Done")
    return {
        "table": table,
        "n_seeds_mapped": n_seeds_mapped,
        "n_seeds_not_in_string": n_not_in_string,
        "n_hits": n_hits,
        "n_recruited": n_recruited,
        "seed_string_ids": seed_string_ids,
        # Diagnostics only, never p-values (see _propagation_core's docstring):
        # how the normal-tail approximation held up against the actual
        # permutation distribution.
        "null_diagnostics": {
            "n_exceed_top": core["n_exceed_top"],
            "null_skew_at_max": core["null_skew_at_max"],
            "n_permutations": core["n_permutations"],
        },
    }


# ---------------------------------------------------------------------------
# Module detection (Markov Clustering / MCL)
# ---------------------------------------------------------------------------

def mcl_cluster(
    adjacency,
    node_ids: list,
    inflation: float = DEFAULT_MCL_INFLATION,
    expansion: int = MCL_EXPANSION,
    max_iter: int = MCL_MAX_ITER,
    prune_threshold: float = MCL_PRUNE_THRESHOLD,
) -> list:
    """Markov Clustering on a small (hit-induced) weighted subnetwork - the
    domain standard for PPI-network module detection (the same family of
    algorithm STRING's own website uses), chosen because it's weight-aware
    (STRING confidence scores directly influence cluster membership) and
    needs no external dependency: it's alternating matrix expansion (spreads
    flow through paths) and inflation (elementwise power + renormalize,
    sharpening strong flows and suppressing weak ones) until convergence.

    adjacency: sparse (n x n), symmetric, nonnegative weights. (A dense array
    is accepted and converted to sparse.)
    node_ids: length-n list identifying each row/col (e.g. STRING protein ids).
    Practical for the SMALL hit subnetwork here (not the whole organism graph).

    Implemented entirely on scipy.sparse (expansion = sparse @ sparse matmul).
    This deliberately avoids dense numpy BLAS matmul (np.linalg.matrix_power /
    dense `@`), which hard-crashes on some Windows numpy/MKL installs when
    Python is launched outside an IDE - sparse matmul is unaffected, and is
    also the standard, more memory-efficient way to run MCL.

    Returns a list of clusters, each a list of node_ids.
    """
    if sp is None:
        raise RuntimeError("scipy is required for MCL clustering")
    n = len(node_ids)
    if n == 0:
        return []
    if n == 1:
        return [[node_ids[0]]]

    M = sp.csr_matrix(adjacency, dtype=np.float64)
    # Self-loops so weakly-connected nodes survive normalization: set each
    # diagonal entry to max(existing, 1).
    diag_add = np.maximum(0.0, 1.0 - M.diagonal())
    if np.any(diag_add > 0):
        M = (M + sp.diags(diag_add)).tocsr()

    def _col_normalize(mat, _sp=sp):
        # _sp=sp (evaluated in the enclosing scope at def-time, where the
        # `if sp is None: raise` guard above has already narrowed it to
        # non-Optional) - a nested function's own body doesn't inherit that
        # narrowing on the closed-over `sp` directly, since Pylance can't
        # prove it wasn't reassigned before the closure is later called.
        sums = np.asarray(mat.sum(axis=0)).ravel()
        sums = np.where(sums == 0, 1.0, sums)
        return _sp.csr_matrix(mat @ _sp.diags(1.0 / sums))

    M = _col_normalize(M)
    for _ in range(max_iter):
        M_prev = M
        # Expansion: matrix power via repeated SPARSE matmul (not dense BLAS).
        expanded = M
        for _e in range(expansion - 1):
            expanded = expanded @ M
        M = expanded.tocsr()
        # Inflation: elementwise power, prune tiny flows, renormalize columns.
        M = M.power(inflation)
        if M.nnz:
            M.data[M.data < prune_threshold] = 0.0
            M.eliminate_zeros()
        M = _col_normalize(M)
        diff = abs(M - M_prev)
        if (diff.max() if diff.nnz else 0.0) < MCL_CONVERGE_TOL:
            break

    # Interpret the converged column-stochastic matrix (dense read-out only -
    # .toarray() copies memory, it is NOT a BLAS matmul, so it's crash-safe).
    # M[i, j] = flow from node j to node i; at convergence node j's flow
    # concentrates on its "attractor" row(s) i. Nodes whose columns flow to a
    # shared attractor belong together, so union each node j with the row(s)
    # its column points to; union-find then merges everything transitively
    # connected (co-attracted nodes → one cluster; a node that only self-
    # attracts stays a singleton).
    Md = M.toarray()
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for j in range(n):
        for i in np.nonzero(Md[:, j] > 0)[0]:
            union(j, int(i))

    groups: dict = {}
    for i in range(n):
        root = find(i)
        groups.setdefault(root, []).append(node_ids[i])
    return list(groups.values())


def detect_modules(
    hit_string_ids: list,
    graph: dict,
    min_module_size: int = DEFAULT_MIN_MODULE_SIZE,
    inflation: float = DEFAULT_MCL_INFLATION,
) -> list:
    """Cluster the hit-induced subnetwork into modules via MCL, dropping
    modules smaller than min_module_size. Returns
    [{"module_id": "M1", "string_ids": [...], "size": int}, ...], largest
    first (a stable, presentation-friendly order)."""
    if sp is None:
        raise RuntimeError("scipy is required for module detection")
    raw_adj = graph["raw_adjacency"]
    protein_index = graph["protein_index"]
    idx = [protein_index[pid] for pid in hit_string_ids if pid in protein_index]
    if not idx:
        return []
    sub = raw_adj[np.ix_(idx, idx)]
    clusters = mcl_cluster(sub, [hit_string_ids[i] for i in range(len(hit_string_ids)) if hit_string_ids[i] in protein_index], inflation=inflation)
    clusters = [c for c in clusters if len(c) >= min_module_size]
    clusters.sort(key=len, reverse=True)
    return [
        {"module_id": f"M{i+1}", "string_ids": members, "size": len(members)}
        for i, members in enumerate(clusters)
    ]


# ---------------------------------------------------------------------------
# Per-module functional enrichment (reuses the ORA Fisher/BH engine)
# ---------------------------------------------------------------------------

def enrich_modules(
    modules: list,
    graph: dict,
    taxon_id: str,
    taxon_name: str = "",
    database: str = "GO",
    include_iea: bool = True,
    min_term_size: int = genomics_ora.MIN_TERM_SIZE,
    max_term_size: int = DEFAULT_MODULE_MAX_TERM_SIZE,
    verbose: bool = False,
) -> dict:
    """Run ORA-style enrichment (Fisher's exact + BH-FDR, reusing
    Genomics_module._run_fisher_for_namespace) on each module's genes, to
    label modules with the functional term that best explains them - the
    step that turns "here are some clusters" into "Module 1 = DNA repair".
    Background/universe is the whole annotated genome (documented; a
    per-module background isn't meaningful here). max_term_size excludes
    huge, uninformative catch-all terms (e.g. "protein binding", ~10k genes)
    from being tested at all for module labeling (via
    _run_fisher_for_namespace's own size gate, same as min_term_size) - without
    it, a module's top term is often such a catch-all instead of a specific,
    interpretable label.

    taxon_name is REQUIRED for database="Reactome" (build_reactome_universe
    filters by organism NAME, unlike GO which filters by taxon_id) - a caller
    that omits it gets every organism's pathways pooled together silently.
    Found live: 91,531 accessions pooled across all species vs 12,155 for
    human alone. GO is unaffected (build_go_universe filters by taxon_id).

    Returns {"results": {module_id: {"namespace": str, "block": <ORA result
    block>} | None}, "caveat": str | None}. The caveat (e.g. a transient
    QuickGO fetch failure, or a GO-hierarchy propagation failure - see
    build_go_universe) applies to every module's enrichment, since they all
    share the same universe, so it is returned once at the top level rather
    than duplicated onto each module's block.
    """
    string_to_uniprot = {v: k for k, v in graph["uniprot_to_protein"].items()}

    if database == "GO":
        universe = genomics_ora.build_go_universe(taxon_id, include_iea=include_iea, verbose=verbose)
        term_maps = {ns: universe[ns] for ns in genomics_ora.GO_NAMESPACES}
    elif database == "Reactome":
        universe = genomics_ora.build_reactome_universe(taxon_id, taxon_name=taxon_name, include_iea=include_iea, verbose=verbose)
        term_maps = {"Reactome": universe["pathways"]}
    else:
        raise ValueError(f"Unsupported module-enrichment database: {database!r}")

    # Precompute each namespace's universe ONCE (not per module) - unioning all
    # term genes is expensive and identical for every module.
    universe_by_ns = {}
    for key, term_map in term_maps.items():
        u: set = set()
        for info in term_map.values():
            u |= info["genes"]
        universe_by_ns[key] = u

    results: dict = {}
    for mod in modules:
        accs = {string_to_uniprot[sid] for sid in mod["string_ids"] if sid in string_to_uniprot}
        if not accs:
            results[mod["module_id"]] = None
            continue
        best_block, best_key = None, None
        for key, term_map in term_maps.items():
            block = genomics_ora._run_fisher_for_namespace(
                term_map, accs, universe_by_ns[key], include_bonferroni=True,
                min_term_size=min_term_size, max_term_size=max_term_size,
            )
            block_min_fdr = block["table"]["FDR"].min() if not block["table"].empty else None
            best_min_fdr = best_block["table"]["FDR"].min() if best_block is not None and not best_block["table"].empty else None
            if best_block is None or (block_min_fdr is not None and (best_min_fdr is None or block_min_fdr < best_min_fdr)):
                best_block, best_key = block, key
        results[mod["module_id"]] = {"namespace": best_key, "block": best_block} if best_block is not None else None
    return {"results": results, "caveat": universe.get("caveat")}
