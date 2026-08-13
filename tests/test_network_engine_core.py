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
print("=== _welford_update: matches numpy's mean/variance ===")
# =============================================================================
rng = np.random.default_rng(7)
samples = [rng.normal(size=5) for _ in range(200)]
mean = np.zeros(5)
m2 = np.zeros(5)
for c, x in enumerate(samples, start=1):
    mean, m2 = gnet._welford_update(mean, m2, c, x)
stack = np.vstack(samples)
check(np.allclose(mean, stack.mean(axis=0), atol=1e-12), "Welford mean matches np.mean to 1e-12")
check(np.allclose(m2 / (len(samples) - 1), stack.var(axis=0, ddof=1), atol=1e-12),
      "Welford m2/(n-1) matches np.var(ddof=1) to 1e-12")

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
    """
    n_nodes = 80
    A = np.zeros((n_nodes, n_nodes))
    W = 0.9

    def clique(nodes):
        for x in range(len(nodes)):
            for y in range(x + 1, len(nodes)):
                A[nodes[x], nodes[y]] = A[nodes[y], nodes[x]] = W

    cluster_a = list(range(0, 10))
    sat_a = [10, 11, 12]
    cluster_b = list(range(13, 23))
    sat_b = [23, 24, 25]
    clique(cluster_a)
    clique(cluster_b)
    for k, s in enumerate(sat_a):                      # each satellite -> 2 cluster members
        for t in (cluster_a[2 * k], cluster_a[2 * k + 1]):
            A[s, t] = A[t, s] = W
    for k, s in enumerate(sat_b):
        for t in (cluster_b[2 * k], cluster_b[2 * k + 1]):
            A[s, t] = A[t, s] = W
    filler = list(range(26, 80))                       # 54 nodes -> 18 disjoint triangles
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
    pg, {f"U{i}": 1.0 for i in seeded_cluster}, n_permutations=300, fdr_alpha=0.05, seed_rng=3)
ptab = pres["table"].set_index("STRING ID")
z_seeded = [float(ptab.loc[f"{TAX}.N{i}", "Z-score"]) for i in sat_seeded]
z_control = [float(ptab.loc[f"{TAX}.N{i}", "Z-score"]) for i in sat_control]
rec = sum(1 for i in sat_seeded if ptab.loc[f"{TAX}.N{i}", "FDR"] < 0.05)
ctrl = sum(1 for i in sat_control if ptab.loc[f"{TAX}.N{i}", "FDR"] < 0.05)

# Separation is the honest measure of power here, not a raw count of FDR<0.05.
# Measured: seeded satellites Z = +3.75/+2.76/+2.70 against controls at
# -1.99/-2.38/-2.15 - unambiguous. Only one clears FDR<0.05 because BH is
# applied across all 80 proteins of which just three are true positives, which
# is a property of the multiple-testing correction at this fixture's size, not
# a weakness of the propagation. Asserting separation keeps the test meaningful
# without hard-coding a threshold that would drift with n.
check(min(z_seeded) > max(z_control),
      f"EVERY satellite of the seeded cluster out-scores EVERY satellite of the unseeded one "
      f"(seeded Z {[round(z,2) for z in z_seeded]} vs control Z {[round(z,2) for z in z_control]}) - "
      f"identical degree and identical local shape, so the signal is topological proximity to seeds, "
      f"not degree")
check(rec >= 1,
      f"at least one seeded satellite is genuinely significant at FDR<0.05 ({rec}/3) - real significance is "
      f"attainable, so the null-calibration result above is not merely 'it never calls anything'")
check(ctrl == 0,
      f"and NO satellite of the unseeded cluster is ever called ({ctrl}/3) - specificity")
check(pres["n_recruited"] >= 1,
      f"the run reports recruited hits, i.e. significant proteins that were never seeds "
      f"(n_recruited={pres['n_recruited']}) - the capability the Transcriptomics module is being built for")

# Nodes no seed can ever reach exercise the null_std == 0 guard: their z-score
# is forced to exactly 0 (conservative) rather than becoming NaN or +inf.
z_filler = [float(ptab.loc[f"{TAX}.N{i}", "Z-score"]) for i in (26, 30, 40, 50)]
check(all(z == 0.0 for z in z_filler),
      f"proteins unreachable from any seed get a z-score of exactly 0, not NaN/inf (got {z_filler}) - "
      f"the null_std==0 guard is real and exercised")

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

TMP.cleanup()
print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
