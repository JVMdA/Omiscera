"""Characterization tests for the SHARED network engine in Genomics_network.py.

WHY THIS FILE EXISTS. Genomics_network.py is a 737-line numerical engine -
STRING parsing, random-walk-with-restart, a degree-preserving permutation
null, MCL clustering - and until this file it had ZERO unit coverage. Nothing
tested _parse_string_network_files, propagate_rwr, the permutation null,
mcl_cluster or detect_modules. The only test touching the module at all
(test_reactome_module_organism_filter.py) stubs the graph out entirely, and
playwright_genomics_results_cards.py asserts nothing about it but a CSS class.

That mattered enough on its own, but it became urgent because the new
Transcriptomics Network / Module Detection module is built directly on these
primitives AND refactors them (extracting an identifier-agnostic propagation
core). Refactoring untested numerical code is how silent wrongness ships. So
these tests are written FIRST, against the CURRENT behaviour, and include a
GOLDEN hit table the refactor must reproduce byte-identically.

Everything here runs on tiny hand-built STRING fixture files in a temporary
directory - no downloads, no server, no network.
"""
import gzip
import os
import sys
import tempfile

import numpy as np
import scipy.sparse as spx
import scipy.sparse.linalg as spla

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import Genomics_network as gnet  # noqa: E402

PASS = FAIL = 0
GOLDEN_PATH = os.path.join(_ROOT, "tests", "fixtures", "golden_netprop_hits.csv")


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


# ---------------------------------------------------------------------------
# Fixtures: a 6-protein graph as two triangles joined by ONE bridge edge whose
# combined_score sits EXACTLY on the default confidence cutoff (700), so the
# boundary rule (>= vs >) is directly observable. P7 appears in the info file
# but in no surviving edge - it must be known by name yet absent from the graph.
# ---------------------------------------------------------------------------
TAX = "9999"
P = {k: f"{TAX}.P{k}" for k in range(1, 8)}

LINKS_ROWS = [
    (P[1], P[2], 900), (P[2], P[1], 900),
    (P[1], P[3], 900), (P[3], P[1], 900),
    (P[2], P[3], 900), (P[3], P[2], 900),
    (P[4], P[5], 900), (P[5], P[4], 900),
    (P[4], P[6], 900), (P[6], P[4], 900),
    (P[5], P[6], 900), (P[6], P[5], 900),
    (P[3], P[4], 700), (P[4], P[3], 700),   # bridge, EXACTLY at the cutoff
    (P[1], P[6], 150), (P[6], P[1], 150),   # below any cutoff we use
]
INFO_ROWS = [(P[i], f"GENE{i}", "100", "annotation") for i in range(1, 8)]
ALIAS_ROWS = [
    (P[1], "Q00001", "UniProt_AC"),
    (P[1], "ENSG00000001", "Ensembl_gene"),      # must NOT enter uniprot_to_protein
    (P[2], "Q00002", "UniProt_AC"),
    (P[3], "Q00003", "UniProt_AC"),
    (P[4], "Q00004", "UniProt_AC"),
    (P[5], "Q00005", "UniProt_AC"),
    (P[6], "Q00006", "UniProt_AC"),
    (P[2], "12345", "UniProt_DR_GeneID"),        # must NOT enter uniprot_to_protein
]


def _write_gz(path, header, rows, sep):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        fh.write(sep.join(header) + "\n")
        for r in rows:
            fh.write(sep.join(str(x) for x in r) + "\n")


def make_fixture_files(tmpdir):
    links = os.path.join(tmpdir, "links.txt.gz")
    info = os.path.join(tmpdir, "info.txt.gz")
    aliases = os.path.join(tmpdir, "aliases.txt.gz")
    _write_gz(links, ["protein1", "protein2", "combined_score"], LINKS_ROWS, " ")
    _write_gz(info, ["#string_protein_id", "preferred_name", "protein_size", "annotation"], INFO_ROWS, "\t")
    _write_gz(aliases, ["#string_protein_id", "alias", "source"], ALIAS_ROWS, "\t")
    return links, info, aliases


TMP = tempfile.TemporaryDirectory()
LINKS, INFO, ALIASES = make_fixture_files(TMP.name)

# =============================================================================
print("=== _parse_string_network_files: structure, boundary, alias filtering ===")
# =============================================================================
g700 = gnet._parse_string_network_files(LINKS, INFO, ALIASES, confidence=700, taxon_id=TAX)
g701 = gnet._parse_string_network_files(LINKS, INFO, ALIASES, confidence=701, taxon_id=TAX)

check(g700["n_proteins"] == 6, f"all 6 edge-bearing proteins are in the graph at confidence 700 (got {g700['n_proteins']})")

i3, i4 = g700["protein_index"][P[3]], g700["protein_index"][P[4]]
bridge_700 = float(g700["raw_adjacency"][i3, i4])
# The pair of real calls below IS the boundary sensitivity check: an edge whose
# score equals the cutoff must be KEPT (>=), and the same edge must vanish one
# unit above it. If the filter were `>`, the first assertion fails immediately.
check(abs(bridge_700 - 0.7) < 1e-12,
      f"an edge scoring EXACTLY the confidence cutoff is KEPT (raw weight {bridge_700}, expected 0.7) - "
      f"the filter is >=, not >")
has_bridge_701 = (P[3] in g701["protein_index"] and P[4] in g701["protein_index"]
                  and float(g701["raw_adjacency"][g701["protein_index"][P[3]], g701["protein_index"][P[4]]]) > 0)
check(not has_bridge_701,
      "and the same edge is GONE one unit above the cutoff - proving the check above is not vacuous, "
      "using two real calls rather than a simulated filter")

diff = abs(g700["raw_adjacency"] - g700["raw_adjacency"].T)
check((diff.max() if diff.nnz else 0.0) < 1e-12, "raw_adjacency is symmetric")

col_sums = np.asarray(g700["adjacency"].sum(axis=0)).ravel()
check(np.allclose(col_sums, 1.0, atol=1e-12),
      f"every column of the normalized adjacency sums to 1 (min {col_sums.min():.6f}, max {col_sums.max():.6f})")

idx_ok = all(g700["index_protein"][i] == pid for pid, i in g700["protein_index"].items())
check(idx_ok, "protein_index and index_protein are mutual inverses")

check(P[7] not in g700["protein_index"] and g700["protein_gene"].get(P[7]) == "GENE7",
      "a protein present in info but in no surviving edge is absent from the graph yet still known by name "
      "(this is what makes a missing 'Gene Symbol' possible downstream)")

check(set(g700["uniprot_to_protein"].keys()) == {"Q00001", "Q00002", "Q00003", "Q00004", "Q00005", "Q00006"},
      f"uniprot_to_protein holds ONLY UniProt_AC rows - Ensembl/GeneID aliases are excluded "
      f"(got {sorted(g700['uniprot_to_protein'].keys())})")

# =============================================================================
print("=== propagate_rwr: conservation, closed form, degenerate input ===")
# =============================================================================
W = g700["adjacency"]
n = g700["n_proteins"]
seed_vec = np.zeros(n)
seed_vec[g700["protein_index"][P[1]]] = 1.0

p = gnet.propagate_rwr(W, seed_vec, restart_prob=0.5)
check(abs(p.sum() - 1.0) < 1e-10, f"the propagated vector is a probability distribution, sums to 1 (got {p.sum():.12f})")

p_alpha1 = gnet.propagate_rwr(W, seed_vec, restart_prob=1.0)
check(np.allclose(p_alpha1, seed_vec / seed_vec.sum(), atol=1e-12),
      "restart_prob=1.0 reproduces the normalized seed exactly (no diffusion at all)")

# Closed form: p = a (I - (1-a) W)^-1 s.
# Solved with scipy.sparse.linalg.spsolve, NOT np.linalg.solve: dense
# numpy/LAPACK calls hard-crash this environment (Windows fatal exception
# 0xc06d007f) when Python runs outside an IDE - the exact failure
# mcl_cluster's own docstring documents and designs around. Sparse solves
# are unaffected, and confirmed working here.
alpha = 0.5
s_norm = seed_vec / seed_vec.sum()
sys_mat = spx.eye(n, format="csc") - (1 - alpha) * spx.csc_matrix(W)
p_closed = alpha * spla.spsolve(sys_mat, s_norm)
check(np.allclose(p, p_closed, atol=1e-10),
      f"the iterative result matches the closed form a(I-(1-a)W)^-1 s to 1e-10 "
      f"(max abs diff {np.abs(p - p_closed).max():.3e})")

i1, i2, i5 = (g700["protein_index"][P[1]], g700["protein_index"][P[2]], g700["protein_index"][P[5]])
check(p[i2] > p[i5],
      f"a direct neighbour of the seed scores above a distant node (P2 {p[i2]:.5f} > P5 {p[i5]:.5f})")

try:
    gnet.propagate_rwr(W, np.zeros(n), restart_prob=0.5)
    check(False, "an all-zero seed vector raises rather than dividing by zero")
except ValueError:
    check(True, "an all-zero seed vector raises ValueError rather than dividing by zero")

# =============================================================================
print("=== _welford_update: matches numpy's mean/variance/skewness ===")
# =============================================================================
rng = np.random.default_rng(7)
# Gamma samples are genuinely right-skewed - unlike a normal sample, this
# actually exercises the M3 term: a bug in the skewness recurrence would not
# necessarily show up on a symmetric distribution, where skew hovers near 0
# regardless of whether M3 is computed correctly.
samples = [rng.gamma(shape=2.0, scale=1.0, size=5) for _ in range(300)]
mean = np.zeros(5)
m2 = np.zeros(5)
m3 = np.zeros(5)
for c, x in enumerate(samples, start=1):
    mean, m2, m3 = gnet._welford_update(mean, m2, m3, c, x)
stack = np.vstack(samples)
check(np.allclose(mean, stack.mean(axis=0), atol=1e-10), "Welford mean matches np.mean to 1e-10")
check(np.allclose(m2 / (len(samples) - 1), stack.var(axis=0, ddof=1), atol=1e-10),
      "Welford m2/(n-1) matches np.var(ddof=1) to 1e-10")
n_s = len(samples)
skew_online = np.sqrt(n_s) * m3 / np.power(m2, 1.5)
sd_pop = stack.std(axis=0, ddof=0)
skew_batch = (((stack - stack.mean(axis=0)) / sd_pop) ** 3).mean(axis=0)
check(np.allclose(skew_online, skew_batch, atol=1e-8),
      f"Welford's streaming skewness (sqrt(n)*M3/M2**1.5) matches a batch computation to 1e-8 "
      f"(max diff {np.abs(skew_online - skew_batch).max():.2e})")
check(skew_batch.mean() > 0.3,
      f"the gamma fixture is genuinely right-skewed (mean skew {skew_batch.mean():.3f} across the 5 "
      f"series), so the check above is not vacuously comparing two near-zero numbers")

# =============================================================================
print("=== _degree_bins: a real partition, ordered by degree ===")
# =============================================================================
bin_of, bins_map = gnet._degree_bins(g700["raw_adjacency"], n_bins=3)
all_binned = sorted(i for members in bins_map.values() for i in members)
check(all_binned == list(range(n)), "every protein index lands in exactly one bin (a true partition)")
degrees = np.asarray(g700["raw_adjacency"].sum(axis=0)).ravel()
bin_mean_deg = [float(np.mean([degrees[i] for i in bins_map[b]])) for b in sorted(bins_map)]
check(all(bin_mean_deg[k] <= bin_mean_deg[k + 1] + 1e-12 for k in range(len(bin_mean_deg) - 1)),
      f"bins are ordered by increasing degree (mean degree per bin: {[round(x,3) for x in bin_mean_deg]})")

# =============================================================================
# Graph builders for the statistical sections. The 6-protein file fixture above
# is deliberately NOT used here: see the degeneracy characterization below.
# =============================================================================
def _graph_from_dense(A, prefix="N"):
    """Wrap a dense symmetric weight matrix in the graph dict shape."""
    n_nodes = A.shape[0]
    raw = spx.csr_matrix(A)
    cs = np.asarray(raw.sum(axis=0)).ravel()
    cs = np.where(cs == 0, 1.0, cs)
    ids = [f"{TAX}.{prefix}{i}" for i in range(n_nodes)]
    return {
        "adjacency": spx.csr_matrix(raw @ spx.diags(1.0 / cs)),
        "raw_adjacency": raw,
        "protein_index": {pid: i for i, pid in enumerate(ids)},
        "index_protein": ids,
        "protein_gene": {pid: f"G{i}" for i, pid in enumerate(ids)},
        "uniprot_to_protein": {f"U{i}": pid for i, pid in enumerate(ids)},
        "n_proteins": n_nodes,
    }


def random_graph(n_nodes, p_edge, rng_):
    A = np.zeros((n_nodes, n_nodes))
    for a in range(n_nodes):
        for b in range(a + 1, n_nodes):
            if rng_.random() < p_edge:
                A[a, b] = A[b, a] = 0.5 + 0.5 * rng_.random()
    return _graph_from_dense(A)


def planted_graph():
    """Two structurally EQUIVALENT high-degree clusters, A and B, each with
    three low-degree satellites. Only A is seeded.

    This shape is required, not incidental. A naive "one dense community, seed
    most of it, recover the rest" fixture cannot show power here, because it
    makes community membership perfectly predict degree - and the null is
    degree-preserving, so every permuted seed set lands back inside the same
    community and the null rises to meet the observed. Giving B the same degree
    profile as A is what lets permuted seeds actually leave the real cluster.

    Expected: A's satellites (adjacent to real seeds) become significant;
    B's satellites (identical degree, identical local shape, but NOT adjacent
    to seeds) must not - that is the specificity half of the test.

    Cluster size (20) and satellite-to-cluster links (4) were sized empirically
    against the SKEW-AWARE p-value (Fix 1), not chosen arbitrarily: the
    original, smaller fixture (10-member clusters, 2 links/satellite) produced
    seeded-satellite Z-scores of ~2.7-3.8, which cleared FDR<0.05 under the
    old (anti-conservative) normal-tail approximation but NOT under the more
    accurate gamma tail - the corrected p-value is legitimately larger for a
    right-skewed null at a given Z, which is the entire point of Fix 1. This
    size was swept (tests/README.txt-adjacent scratch work, not committed)
    across several (cluster_size, links, seed_rng) combinations and gives a
    comfortable, reproducible margin (FDR <= 0.007 for every seeded satellite
    across 3 different seed_rng values), not a value sitting right at the
    0.05 boundary.
    """
    n_nodes = 100
    A = np.zeros((n_nodes, n_nodes))
    W = 0.9

    def clique(nodes):
        for x in range(len(nodes)):
            for y in range(x + 1, len(nodes)):
                A[nodes[x], nodes[y]] = A[nodes[y], nodes[x]] = W

    cluster_a = list(range(0, 20))
    sat_a = [20, 21, 22]
    cluster_b = list(range(23, 43))
    sat_b = [43, 44, 45]
    clique(cluster_a)
    clique(cluster_b)
    for k, s in enumerate(sat_a):                      # each satellite -> 4 cluster members
        for j in range(4):
            t = cluster_a[(2 * k + j) % len(cluster_a)]
            A[s, t] = A[t, s] = W
    for k, s in enumerate(sat_b):
        for j in range(4):
            t = cluster_b[(2 * k + j) % len(cluster_b)]
            A[s, t] = A[t, s] = W
    filler = list(range(46, 100))                      # 54 nodes -> 18 disjoint triangles
    for g in range(0, len(filler), 3):
        clique(filler[g:g + 3])
    return _graph_from_dense(A), cluster_a, sat_a, sat_b


# =============================================================================
print("=== the permutation null: determinism and DEGREE PRESERVATION ===")
# =============================================================================
perm_graph = random_graph(60, 0.12, np.random.default_rng(99))
perm_seeds = {f"U{i}": 1.0 for i in (3, 7, 11, 19, 28, 41)}

r_a = gnet.run_propagation_with_significance(perm_graph, perm_seeds, n_permutations=25, seed_rng=42)
r_b = gnet.run_propagation_with_significance(perm_graph, perm_seeds, n_permutations=25, seed_rng=42)
check(r_a["table"].equals(r_b["table"]), "the same seed_rng reproduces an identical result table (deterministic)")
r_c = gnet.run_propagation_with_significance(perm_graph, perm_seeds, n_permutations=25, seed_rng=1234)
check(not r_a["table"]["Z-score"].equals(r_c["table"]["Z-score"]),
      "a different seed_rng genuinely changes the null (so seed_rng is actually wired through)")

# Degree preservation: intercept the REAL propagate_rwr to capture every
# permuted seed vector, then compare the degree-bin multiset of permuted seed
# positions against the observed seeds'. This is the load-bearing property -
# without it a hub would look significant merely for being a hub.
_src_before = open(gnet.__file__, "rb").read()
captured = []
_real_rwr = gnet.propagate_rwr


def _spy_rwr(adjacency, seed_vector, **kw):
    captured.append(np.nonzero(np.asarray(seed_vector) > 0)[0].tolist())
    return _real_rwr(adjacency, seed_vector, **kw)


gnet.propagate_rwr = _spy_rwr
try:
    gnet.run_propagation_with_significance(perm_graph, perm_seeds, n_permutations=60, seed_rng=5)
finally:
    gnet.propagate_rwr = _real_rwr
check(open(gnet.__file__, "rb").read() == _src_before,
      "the module's source file on disk was never touched by the spy (in-memory monkeypatch only)")

bin_of_full, _ = gnet._degree_bins(perm_graph["raw_adjacency"])
observed_bins = sorted(int(bin_of_full[i]) for i in captured[0])
perm_bin_sets = [sorted(int(bin_of_full[i]) for i in c) for c in captured[1:]]
check(len(perm_bin_sets) == 60, f"one permutation vector was captured per permutation (got {len(perm_bin_sets)})")
check(all(b == observed_bins for b in perm_bin_sets),
      f"EVERY permuted seed set has exactly the observed seeds' degree-bin composition {observed_bins} "
      f"- the null is degree-preserving, not a uniform reshuffle")
moved = sum(1 for c in captured[1:] if sorted(c) != sorted(captured[0]))
check(moved >= 55,
      f"and the permutations genuinely move the seeds ({moved}/60 differ from the observed set) - "
      f"proving the check above is not vacuously comparing identical sets")

print("=== CHARACTERIZED LIMITATION: the null degenerates on a tiny graph ===")
# Found while writing these tests, and real: _degree_bins uses
# bin_size = max(1, n // DEGREE_BIN_COUNT), so for n < 20 every bin holds at
# most ONE protein. The permutation can then only ever re-pick each seed
# itself, making the null identical to the observed, every z-score 0, and
# nothing significant - silently, with no warning anywhere. It fails SAFE (no
# false positives) but returns a uniformly uninformative result that looks
# like a normal one. Real STRING organism graphs have thousands of proteins,
# but a sparse organism at high confidence could land here. Asserted so the
# behaviour is documented and cannot change unnoticed; the Transcriptomics
# module guards against it explicitly rather than inheriting it.
tiny = random_graph(12, 0.45, np.random.default_rng(4))
captured.clear()
gnet.propagate_rwr = _spy_rwr
try:
    tiny_res = gnet.run_propagation_with_significance(tiny, {"U1": 1.0, "U5": 1.0}, n_permutations=15, seed_rng=8)
finally:
    gnet.propagate_rwr = _real_rwr
tiny_moved = sum(1 for c in captured[1:] if sorted(c) != sorted(captured[0]))
check(tiny_moved == 0,
      f"on a 12-protein graph NO permutation can move a seed ({tiny_moved}/15 differ) - the degree bins hold "
      f"one candidate each")
check(int((tiny_res['table']['FDR'] < 0.05).sum()) == 0,
      "and consequently nothing is ever significant on such a graph - documented here so a future change "
      "to DEGREE_BIN_COUNT or the binning cannot silently alter it")

# =============================================================================
print("=== FDR calibration under a TRUE null (does this engine invent hits?) ===")
# =============================================================================
null_rates = []
for rep in range(5):
    rg = np.random.default_rng(100 + rep)
    ng = random_graph(60, 0.12, rg)
    seed_accs = {f"U{i}": 1.0 for i in rg.choice(60, size=8, replace=False)}
    res = gnet.run_propagation_with_significance(ng, seed_accs, n_permutations=200,
                                                 fdr_alpha=0.05, seed_rng=int(rg.integers(1e6)))
    null_rates.append(float((res["table"]["FDR"] < 0.05).mean()))
worst = max(null_rates)
check(worst <= 0.05,
      f"under a TRUE null (random graph, random seeds, no planted signal) the FDR<0.05 call rate stays at or "
      f"below 5% across 5 repeats (worst {worst:.3f}, rates {[round(r,3) for r in null_rates]}) - "
      f"the engine does not manufacture significance")

# =============================================================================
print("=== power AND specificity on a planted structure ===")
# =============================================================================
pg, seeded_cluster, sat_seeded, sat_control = planted_graph()
pres = gnet.run_propagation_with_significance(
    pg, {f"U{i}": 1.0 for i in seeded_cluster}, n_permutations=500, fdr_alpha=0.05, seed_rng=3)
ptab = pres["table"].set_index("STRING ID")
z_seeded = [float(ptab.loc[f"{TAX}.N{i}", "Z-score"]) for i in sat_seeded]
z_control = [float(ptab.loc[f"{TAX}.N{i}", "Z-score"]) for i in sat_control]
rec = sum(1 for i in sat_seeded if ptab.loc[f"{TAX}.N{i}", "FDR"] < 0.05)
ctrl = sum(1 for i in sat_control if ptab.loc[f"{TAX}.N{i}", "FDR"] < 0.05)

# Separation is the honest measure of power here, not a raw count of FDR<0.05.
# Measured (seed_rng=3, skew-aware p-values): seeded satellites Z =
# +4.05/+4.25/+3.82, FDR = 0.0049/0.0049/0.0068 against controls at
# -3.46/-3.72/-3.63 - unambiguous, and every seeded satellite clears FDR<0.05
# with a comfortable margin (checked, not assumed, across seed_rng in
# {3, 7, 11} while sizing this fixture - see planted_graph()'s docstring).
check(min(z_seeded) > max(z_control),
      f"EVERY satellite of the seeded cluster out-scores EVERY satellite of the unseeded one "
      f"(seeded Z {[round(z,2) for z in z_seeded]} vs control Z {[round(z,2) for z in z_control]}) - "
      f"identical degree and identical local shape, so the signal is topological proximity to seeds, "
      f"not degree")
check(rec == 3,
      f"ALL THREE seeded satellites are genuinely significant at FDR<0.05 ({rec}/3) - real significance is "
      f"attainable under the skew-aware p-value too, so the null-calibration result above is not merely "
      f"'it never calls anything'")
check(ctrl == 0,
      f"and NO satellite of the unseeded cluster is ever called ({ctrl}/3) - specificity")
check(pres["n_recruited"] >= 1,
      f"the run reports recruited hits, i.e. significant proteins that were never seeds "
      f"(n_recruited={pres['n_recruited']}) - the capability the Transcriptomics module is being built for")

# Nodes no seed can ever reach exercise the null_std == 0 guard: their z-score
# is forced to exactly 0 (conservative) rather than becoming NaN or +inf.
z_filler = [float(ptab.loc[f"{TAX}.N{i}", "Z-score"]) for i in (46, 55, 70, 90)]
check(all(z == 0.0 for z in z_filler),
      f"proteins unreachable from any seed get a z-score of exactly 0, not NaN/inf (got {z_filler}) - "
      f"the null_std==0 guard is real and exercised")

# =============================================================================
print("=== Fix 1: the skew-aware p-value is genuinely closer to the empirical tail than a normal-only one ===")
# =============================================================================
# Direct call to _propagation_core (not the identifier-vocabulary wrapper) to
# reach the per-protein diagnostics (null_skew, pvalues_empirical) the wrapper
# does not expose. Reuses the same planted graph and seeds as the power/
# specificity check above, so this is real data, not a synthetic distribution.
seed_idx_score_pg, _ = gnet._map_seeds_to_indices(
    pg["uniprot_to_protein"], pg["protein_index"], {f"U{i}": 1.0 for i in seeded_cluster})
core_pg = gnet._propagation_core(pg, seed_idx_score_pg, n_permutations=500, seed_rng=3)
skew_pg = core_pg["null_skew"]
high_skew_idx = int(np.argmax(skew_pg))
check(skew_pg[high_skew_idx] > gnet.SKEW_GAMMA_THRESHOLD,
      f"the fixture genuinely produces a protein with null skew above SKEW_GAMMA_THRESHOLD "
      f"(max skew {skew_pg[high_skew_idx]:.2f}) - otherwise the gamma fit is never exercised and the "
      f"comparison below would be vacuous")
check(core_pg["n_gamma_tail"] > 0,
      f"at least one protein's p-value actually used the gamma fit, not just the normal fallback "
      f"(n_gamma_tail={core_pg['n_gamma_tail']} of {pg['n_proteins']})")

from scipy.special import ndtr as _ndtr

p_empirical = float(core_pg["pvalues_empirical"][high_skew_idx])
p_skew_aware = float(core_pg["pvalues"][high_skew_idx])
z_here = float(core_pg["zscores"][high_skew_idx])
p_normal_only = float(_ndtr(-z_here))
log_err_skew = abs(np.log10(max(p_skew_aware, 1e-300)) - np.log10(max(p_empirical, 1e-300)))
log_err_normal = abs(np.log10(max(p_normal_only, 1e-300)) - np.log10(max(p_empirical, 1e-300)))
check(log_err_skew < log_err_normal,
      f"for the highest-skew protein (skew {skew_pg[high_skew_idx]:.2f}), the skew-aware p-value "
      f"({p_skew_aware:.3g}) sits closer in log-space to the empirical exceedance rate "
      f"({p_empirical:.3g}) than the pre-fix normal-only p-value ({p_normal_only:.3g}) does "
      f"(log10 error {log_err_skew:.2f} vs {log_err_normal:.2f}) - Fix 1 measurably improves calibration "
      f"here, not just in theory")

# =============================================================================
print("=== Fix 2: fdrs_recruited_family is a real, independent BH family over non-seed p-values ===")
# =============================================================================
seed_positions = np.fromiter(seed_idx_score_pg.keys(), dtype=int)
non_seed_positions = np.setdiff1d(np.arange(pg["n_proteins"]), seed_positions)
fdrs_rf = core_pg["fdrs_recruited_family"]
check(np.all(np.isnan(fdrs_rf[seed_positions])),
      f"every SEED position is NaN in fdrs_recruited_family - it is not a question for a seed "
      f"(got {np.sum(~np.isnan(fdrs_rf[seed_positions]))} non-NaN seed entries, expected 0)")
check(not np.any(np.isnan(fdrs_rf[non_seed_positions])),
      "every NON-seed position has a real value")

manual_bh = np.asarray(gnet.genomics_ora.benjamini_hochberg(
    [float(core_pg["pvalues"][i]) for i in non_seed_positions]))
check(np.allclose(fdrs_rf[non_seed_positions], manual_bh, atol=1e-12),
      "fdrs_recruited_family matches BH benjamini_hochberg recomputed independently over just the "
      "SAME non-seed p-values (max abs diff %.2e)" % np.max(np.abs(fdrs_rf[non_seed_positions] - manual_bh)))

# SENSITIVITY: prove the recruited-only family is not merely an alias for the
# all-protein "fdrs" column read for non-seed rows - it must differ, since
# removing every seed's (structurally tiny) p-value from the correction
# changes the family size AND the rank of every remaining candidate.
all_family_for_non_seed = np.asarray(core_pg["fdrs"])[non_seed_positions]
check(not np.allclose(fdrs_rf[non_seed_positions], all_family_for_non_seed, atol=1e-9),
      f"the recruited-only family FDR genuinely differs from the all-protein family FDR for the SAME "
      f"non-seed proteins (max abs diff {np.max(np.abs(fdrs_rf[non_seed_positions] - all_family_for_non_seed)):.3e}) "
      f"- confirming this is a real second computation, not a relabelled copy of the first")

print("=== SENSITIVITY: the fdrs_recruited_family check above is not vacuous ===")
# Break it for real: compute the "recruited family" from ALL p-values (the
# bug this fix replaces - a naive read of the shared "fdrs" column) instead
# of the non-seed subset, and confirm the match-to-manual-BH check fails.
broken_fdrs_rf = np.asarray(gnet.genomics_ora.benjamini_hochberg(list(core_pg["pvalues"])))[non_seed_positions]
check(not np.allclose(broken_fdrs_rf, manual_bh, atol=1e-12),
      f"the PRE-FIX approach (BH over ALL p-values, not just non-seed) does NOT match the independent "
      f"non-seed-only BH recomputation (max abs diff "
      f"{np.max(np.abs(broken_fdrs_rf - manual_bh)):.3e}) - proving the match check above is discriminating, "
      f"not accidentally always true")

# =============================================================================
print("=== mcl_cluster / detect_modules ===")
# =============================================================================
# Back to the small hand-built file fixture: module detection is pure graph
# structure, so the 6-protein two-triangle graph is the clearest possible case
# (the statistical sections above needed a larger graph only because the
# degree-binned permutation is meaningless at that size - see the degeneracy
# characterization).
graph = g700
two_cliques = np.zeros((6, 6))
for a, b in [(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5)]:
    two_cliques[a, b] = two_cliques[b, a] = 0.9
ids6 = ["A", "B", "C", "D", "E", "F"]
clusters = gnet.mcl_cluster(spx.csr_matrix(two_cliques), ids6, inflation=2.0)
sizes = sorted(len(c) for c in clusters)
check(sizes == [3, 3], f"two disjoint triangles cluster into exactly two 3-member modules (got sizes {sizes})")
check(all(set(c) in ({"A", "B", "C"}, {"D", "E", "F"}) for c in clusters),
      f"and the membership is correct, returning the CALLER's node_ids not indices (got {clusters})")
check(gnet.mcl_cluster(spx.csr_matrix(np.zeros((0, 0))), [], inflation=2.0) == [], "n=0 returns no clusters")
check(gnet.mcl_cluster(spx.csr_matrix(np.zeros((1, 1))), ["solo"], inflation=2.0) == [["solo"]], "n=1 returns one singleton")

# detect_modules: an unknown id in the MIDDLE must not shift id<->index alignment.
mods = gnet.detect_modules([P[1], "9999.NOT_IN_GRAPH", P[2], P[3]], graph, min_module_size=2, inflation=2.0)
member_ids = {m for mod in mods for m in mod["string_ids"]}
check("9999.NOT_IN_GRAPH" not in member_ids, "an unknown STRING id never appears in a module")
check(member_ids <= {P[1], P[2], P[3]} and len(member_ids) == 3,
      f"the three real hits are clustered with correct identities despite the unknown id sitting BETWEEN them "
      f"(got {sorted(member_ids)}) - the two paired comprehensions stay aligned")
check(all(mods[k]["size"] >= mods[k + 1]["size"] for k in range(len(mods) - 1)),
      "modules are returned largest-first")
check([m["module_id"] for m in mods] == [f"M{k+1}" for k in range(len(mods))],
      f"module ids are contiguous M1..Mk (got {[m['module_id'] for m in mods]})")
check(gnet.detect_modules([P[1], P[2], P[3]], graph, min_module_size=99) == [],
      "min_module_size genuinely filters modules out")

# =============================================================================
print("=== GOLDEN: the exact genomics hit table, for the upcoming refactor ===")
# =============================================================================
golden_run = gnet.run_propagation_with_significance(graph, {"Q00001": 2.5, "Q00004": 1.0},
                                                    restart_prob=0.5, n_permutations=50,
                                                    fdr_alpha=0.05, seed_rng=2024)
golden_csv = golden_run["table"].to_csv(index=False)
os.makedirs(os.path.dirname(GOLDEN_PATH), exist_ok=True)
if not os.path.exists(GOLDEN_PATH):
    with open(GOLDEN_PATH, "w", encoding="utf-8", newline="") as fh:
        fh.write(golden_csv)
    print(f"    (golden baseline CREATED at {os.path.relpath(GOLDEN_PATH, _ROOT)} - "
          f"commit it; from now on it is the refactor contract)")
    check(True, "golden baseline written from current pre-refactor behaviour")
else:
    with open(GOLDEN_PATH, "r", encoding="utf-8", newline="") as fh:
        stored = fh.read()
    check(stored == golden_csv,
          "the genomics propagation output is byte-identical to the committed golden baseline - "
          "any refactor of the shared engine has NOT changed genomics behaviour")

check(list(golden_run["table"].columns) == gnet.NETPROP_HIT_COLUMNS,
      f"the hit table's columns are exactly NETPROP_HIT_COLUMNS (got {list(golden_run['table'].columns)})")
check(golden_run["n_seeds_mapped"] == 2 and golden_run["n_seeds_not_in_string"] == 0,
      f"seed accounting is correct (mapped {golden_run['n_seeds_mapped']}, "
      f"unmapped {golden_run['n_seeds_not_in_string']})")
unmapped = gnet.run_propagation_with_significance(graph, {"Q00001": 1.0, "NOT_A_REAL_ACC": 1.0},
                                                   n_permutations=10, seed_rng=1)
check(unmapped["n_seeds_not_in_string"] == 1,
      f"an accession absent from STRING is counted as unmapped, not silently dropped "
      f"(got {unmapped['n_seeds_not_in_string']})")
try:
    gnet.run_propagation_with_significance(graph, {"NOPE": 1.0}, n_permutations=5)
    check(False, "a run where NOTHING maps raises rather than returning noise")
except ValueError:
    check(True, "a run where NOTHING maps raises ValueError rather than returning noise")

# =============================================================================
print("=== _GRAPH_CACHE is thread-safe (API.py's EXECUTOR has 4 worker threads) ===")
# =============================================================================
import threading
import time as _time

# 1. The SAME lock object actually serializes _graph_cache_put calls, not just
# "a lock exists somewhere in the file". A background thread holds the real
# lock for 300ms; a concurrent put() must block until it releases.
hold_started = threading.Event()
release_at = [None]


def _hold_lock():
    with gnet._GRAPH_CACHE_LOCK:
        hold_started.set()
        _time.sleep(0.3)
        release_at[0] = _time.monotonic()


holder = threading.Thread(target=_hold_lock)
holder.start()
hold_started.wait(timeout=2)
t0 = _time.monotonic()
gnet._graph_cache_put(("lock-race-test", 1), {"n_proteins": 0})
put_finished_at = _time.monotonic()
holder.join()
check(release_at[0] is not None and put_finished_at >= release_at[0],
      f"_graph_cache_put genuinely blocked until the holder released the SAME lock "
      f"(put finished at {put_finished_at - t0:.3f}s, holder released at "
      f"{release_at[0] - t0:.3f}s)")
gnet._GRAPH_CACHE.pop(("lock-race-test", 1), None)

# 2. Many threads hammering build_string_network across more distinct keys
# than the cache holds, with the real STRING download/parse stubbed out so
# this stays fast and network-free. A pre-fix OrderedDict here would
# occasionally KeyError on move_to_end racing another thread's popitem.
_orig_ensure_file = gnet._ensure_string_file
_orig_parse = gnet._parse_string_network_files
_call_errors = []
_call_lock = threading.Lock()


def _stub_ensure_file(url, local_path, refresh_days=None, verbose=False):
    return local_path  # never touches disk or network


def _stub_parse(links_path, info_path, aliases_path, confidence, taxon_id=""):
    _time.sleep(0.005)  # widen the window for a real interleaving
    return {"n_proteins": int(taxon_id) if taxon_id.isdigit() else -1, "tag": taxon_id}


gnet._ensure_string_file = _stub_ensure_file
gnet._parse_string_network_files = _stub_parse
try:
    gnet._GRAPH_CACHE.clear()
    taxa = [str(100 + (i % 5)) for i in range(40)]  # 5 distinct keys, cache holds 2

    def _worker(tax):
        try:
            g = gnet.build_string_network(tax, confidence=700, verbose=False)
            if g.get("n_proteins") != int(tax):
                with _call_lock:
                    _call_errors.append(f"wrong graph returned for taxon {tax}: {g}")
        except Exception as exc:
            with _call_lock:
                _call_errors.append(f"{tax}: {exc!r}")

    threads = [threading.Thread(target=_worker, args=(t,)) for t in taxa]
    for th in threads:
        th.start()
    for th in threads:
        th.join(timeout=10)
finally:
    gnet._ensure_string_file = _orig_ensure_file
    gnet._parse_string_network_files = _orig_parse

check(_call_errors == [], f"40 concurrent build_string_network calls across 5 keys, 0 errors "
                          f"(got {len(_call_errors)}: {_call_errors[:3]})")
check(len(gnet._GRAPH_CACHE) <= gnet._GRAPH_CACHE_MAX_ENTRIES,
      f"the cache never exceeded its {gnet._GRAPH_CACHE_MAX_ENTRIES}-entry bound under concurrent "
      f"access (got {len(gnet._GRAPH_CACHE)})")
gnet._GRAPH_CACHE.clear()

TMP.cleanup()
print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
