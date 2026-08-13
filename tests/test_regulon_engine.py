"""Direct unit tests for the new Transcriptomics_regulon.py engine (TF /
Regulon module). No server, no network - pure functions + hand-built
fixtures, mirroring test_burden_engine.py / test_netprop_engine.py's approach.

Priorities, in the order the module docstring lists its own risk:
  1. The ULM closed-form (_ulm_score) against a naive, brute-force OLS fit on
     the SAME fixture, by hand-computed numbers - this is the one genuinely
     new statistic in the whole module.
  2. SIGN correctness end-to-end through run_regulon_activity: a TF whose
     activated targets are all up and repressed targets all down must score
     clearly positive, and flipping the DE table must flip the sign. This is
     exactly the failure mode a mis-read OmniPath direction field would
     produce, and it would look plausible rather than obviously broken.
  3. Degenerate inputs (constant y, zero-variance x, undersized regulon, a TF
     with none of its targets present in the table).
  4. Parsers: TRRUST TSV, ChEA3 GMT (incl. the cell-line-union behaviour),
     OmniPath JSON rows (incl. the COMPLEX: exclusion).
  5. Cache round-trip: sets/weights survive the JSON (de)serialization.
"""
import math
import os
import shutil
import sys
import tempfile

import numpy as np
import pandas as pd

# Repo-relative: this suite lives in <repo>/tests/ and the modules under
# test are one level up. Was a hardcoded absolute path while these were
# throwaway scratchpad scripts; made portable when they moved into the repo.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, _ROOT)
import Transcriptomics_regulon as tr

# Any test below that fetches through a FAKE HTTP session and lets a fetcher
# reach its own force_refresh=True save path (fetch_trrust_regulons /
# fetch_chea3_regulons) would otherwise write fixture data straight into the
# real on-disk cache at tr.REGULON_CACHE_DIR - the SAME file the live app
# reads (cache key is only the taxon/library name, not "is this a test").
# Caught in this session: an earlier run of this exact file wrote fake
# TF1/TF2 TRRUST data over the real trrust_9606.json cache, and the app
# silently served it for the rest of that session (30-day refresh window).
# Every fetcher call below is redirected to a throwaway temp directory for
# the whole file, restored at the end - never assume a "test" fetch is safe
# just because the session is mocked.
_REAL_REGULON_CACHE_DIR = tr.REGULON_CACHE_DIR
_TEST_CACHE_DIR = tempfile.mkdtemp(prefix="omiscera_regulon_test_cache_")
tr.REGULON_CACHE_DIR = _TEST_CACHE_DIR

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


# =============================================================================
# 1. _ulm_score vs. brute-force OLS on a hand-checkable fixture
# =============================================================================
print("=== _ulm_score: closed-form vs. brute-force OLS ===")

# 10 genes; y = ranking metric. Genes 0-2 are "activated targets" (x=+1),
# genes 3-4 are "repressed targets" (x=-1), genes 5-9 are non-targets (x=0).
y = np.array([5.0, 4.0, 6.0, -3.0, -2.0, 0.5, -0.5, 1.0, -1.0, 0.2])
n = len(y)
pos_positions = [0, 1, 2]
neg_positions = [3, 4]

y_mean = y.mean()
yc = y - y_mean
syy = float(np.dot(yc, yc))

score, p_value, n_used, beta1 = tr._ulm_score(pos_positions, neg_positions, yc, syy, n)
check(n_used == 5, f"n_used == 5 (got {n_used})")
check(score is not None, "score computed (non-degenerate fixture)")
check(beta1 is not None, "beta1 returned alongside the score")

# Brute-force reference: build x explicitly, fit OLS by hand (numpy polyfit-
# equivalent normal equations), compute the residual SE and t-statistic
# exactly as a textbook simple-linear-regression t-test would.
x = np.zeros(n)
for i in pos_positions:
    x[i] = 1.0
for i in neg_positions:
    x[i] = -1.0
x_mean = x.mean()
xc = x - x_mean
sxx_ref = float(np.dot(xc, xc))
sxy_ref = float(np.dot(xc, yc))
beta1_ref = sxy_ref / sxx_ref
beta0_ref = y_mean - beta1_ref * x_mean
y_hat = beta0_ref + beta1_ref * x
resid = y - y_hat
sse_ref = float(np.dot(resid, resid))
df_ref = n - 2
se_ref = math.sqrt(sse_ref / df_ref / sxx_ref)
score_ref = beta1_ref / se_ref
from scipy.special import stdtr as _stdtr
p_ref = 2.0 * float(_stdtr(df_ref, -abs(score_ref)))

check(abs(score - score_ref) < 1e-9, f"score matches brute-force OLS t-stat (got {score}, ref {score_ref})")
check(abs(p_value - p_ref) < 1e-9, f"p-value matches brute-force OLS (got {p_value}, ref {p_ref})")
check(abs(beta1 - beta1_ref) < 1e-9, f"beta1 matches brute-force OLS slope (got {beta1}, ref {beta1_ref})")

# Independently hand-verify beta1_ref and sxx_ref too, not just cross-check
# two implementations that could share the same mistake.
# sum(x) = 3*1 + 2*(-1) = 1 -> x_mean = 1/10 = 0.1
check(abs(x_mean - 0.1) < 1e-12, f"x_mean == 0.1 (hand-computed) (got {x_mean})")
# sxx = sum(x^2) - n*x_mean^2 = 5 - 10*0.01 = 4.9
check(abs(sxx_ref - 4.9) < 1e-9, f"sxx == 4.9 (hand-computed) (got {sxx_ref})")


# =============================================================================
# 2. Degenerate inputs
# =============================================================================
print("=== _ulm_score: degenerate inputs ===")

# Fewer than 2 usable targets.
s, p, n_u, b1 = tr._ulm_score([0], [], yc, syy, n)
check(s is None and p is None and b1 is None, "single-target regulon (n_total=1) returns (None, None, ., None)")
check(n_u == 1, f"n_used still reported as 1 (got {n_u})")

s, p, n_u, b1 = tr._ulm_score([], [], yc, syy, n)
check(s is None and n_u == 0, "empty regulon returns (None, None, 0, None)")

# Zero-variance x: equal counts of + and - is fine (sxx>0); the real
# zero-variance case is when EVERY tested gene is a target with the SAME
# sign, which can't happen from real data (target count <= n always), so
# instead verify the sxx<=0 guard directly by constructing an n=2 case where
# x_mean would drive sxx to ~0 - not reachable from real pos/neg counts, so
# this exercises the guard rather than a realistic input.
s, p, n_u, b1 = tr._ulm_score([0, 1], [], np.zeros(2), 0.0, 2)
check(s is None, "zero-target-count-vs-n edge case handled without raising")

# Perfect fit (sse == 0): construct y that is an EXACT linear function of x
# for a tiny n, driving residuals to zero.
y_perfect = np.array([1.0, 1.0, -1.0, -1.0])
yc_perfect = y_perfect - y_perfect.mean()
syy_perfect = float(np.dot(yc_perfect, yc_perfect))
s, p, n_u, b1 = tr._ulm_score([0, 1], [2, 3], yc_perfect, syy_perfect, 4)
check(s is None, "perfect-fit (sse=0) degenerate case returns None rather than raising/inf")


# =============================================================================
# 3. run_regulon_activity: end-to-end sign correctness
# =============================================================================
print("=== run_regulon_activity: sign correctness end-to-end ===")


def _make_resolved_df(log2fc_by_gene, pvalue=0.01):
    rows = []
    for i, (sym, lfc) in enumerate(log2fc_by_gene.items()):
        rows.append({
            "gene_symbol": sym, "entrez_id": str(1000 + i), "ensembl_id": f"ENSG{i:05d}",
            "log2FC": lfc, "pvalue": pvalue, "padj": pvalue,
        })
    return pd.DataFrame(rows)


# TF_UP: activates A1/A2/A3 (should go up), represses R1/R2 (should go down).
# 15 filler genes with near-zero signal complete the "rest of the experiment".
genes_up_scenario = {
    "A1": 3.0, "A2": 2.5, "A3": 2.8,     # activated targets: up
    "R1": -2.5, "R2": -3.0,               # repressed targets: down
}
for i in range(15):
    genes_up_scenario[f"FILLER{i}"] = 0.05 * (-1 if i % 2 else 1)

resolved_up = _make_resolved_df(genes_up_scenario)

fake_regulons = {
    "regulons": {
        "TF1": {
            "name": "TF1",
            "genes": {"A1", "A2", "A3", "R1", "R2"},
            "weights": {"A1": 1, "A2": 1, "A3": 1, "R1": -1, "R2": -1},
        },
    },
    "signed": True, "fetched_at": "2026-01-01T00:00:00+00:00", "fetch_failed": False, "note": None,
}


def _patch_fetch(fake):
    original = tr.fetch_regulons
    tr.fetch_regulons = lambda *a, **k: fake
    return original


original_fetch = _patch_fetch(fake_regulons)
try:
    result_up = tr.run_regulon_activity(
        resolved_up, ["TRRUST"], "9606", metric="log2fc", min_regulon_size=3,
    )
finally:
    tr.fetch_regulons = original_fetch

table_up = result_up["tables"]["TRRUST"]["table"]
check(len(table_up) == 1, f"TF1 present in activity table (got {len(table_up)} rows)")
if len(table_up) == 1:
    score_up = table_up.iloc[0]["Activity Score"]
    check(score_up > 0, f"TF1 scores clearly POSITIVE when its activated targets are up (got {score_up})")

# Flip the DE table (negate every log2FC) - the same TF must now score negative.
resolved_down = _make_resolved_df({g: -v for g, v in genes_up_scenario.items()})
original_fetch = _patch_fetch(fake_regulons)
try:
    result_down = tr.run_regulon_activity(
        resolved_down, ["TRRUST"], "9606", metric="log2fc", min_regulon_size=3,
    )
finally:
    tr.fetch_regulons = original_fetch

table_down = result_down["tables"]["TRRUST"]["table"]
check(len(table_down) == 1, f"TF1 present after flip (got {len(table_down)} rows)")
if len(table_down) == 1:
    score_down = table_down.iloc[0]["Activity Score"]
    check(score_down < 0, f"TF1 scores clearly NEGATIVE after flipping the DE table (got {score_down})")
    check(abs(score_up + score_down) < 1e-6, f"flip is an exact sign flip (up={score_up}, down={score_down})")

# Deliberately break the sign convention (swap + and -) and confirm the
# up-scenario test would then go RED - i.e. the test is actually sensitive to
# this bug class, not just passing by construction.
broken_regulons = {
    "regulons": {
        "TF1": {
            "name": "TF1", "genes": {"A1", "A2", "A3", "R1", "R2"},
            "weights": {"A1": -1, "A2": -1, "A3": -1, "R1": 1, "R2": 1},  # DELIBERATELY BROKEN: signs swapped
        },
    },
    "signed": True, "fetched_at": "x", "fetch_failed": False, "note": None,
}
original_fetch = _patch_fetch(broken_regulons)
try:
    result_broken = tr.run_regulon_activity(resolved_up, ["TRRUST"], "9606", metric="log2fc", min_regulon_size=3)
finally:
    tr.fetch_regulons = original_fetch
broken_score = result_broken["tables"]["TRRUST"]["table"].iloc[0]["Activity Score"]
check(broken_score < 0, f"sanity: swapped-sign fixture DOES score negative (proves the real test is sensitive) (got {broken_score})")


# =============================================================================
# 4. ChEA3 (unsigned) is excluded from activity, present in enrichment path
# =============================================================================
print("=== unsigned database handling ===")
result_mixed = tr.run_regulon_activity(resolved_up, ["ChEA3"], "9606", metric="log2fc")
chea3_block = result_mixed["tables"]["ChEA3"]
check(chea3_block["n_tfs_tested"] == 0, "ChEA3 gets zero TFs tested in activity inference")
check(chea3_block["note"] is not None and "mode-of-regulation" in chea3_block["note"],
      "ChEA3 activity block explains WHY (unsigned) rather than silently being empty")


# =============================================================================
# 5. Parsers
# =============================================================================
print("=== TRRUST TSV parsing ===")
trrust_text = (
    "TF1\tG1\tActivation\t111\n"
    "TF1\tG2\tRepression\t222\n"
    "TF1\tG3\tUnknown\t333\n"
    "TF2\tG1\tActivation\t444\n"
)


class _FakeResp:
    def __init__(self, text):
        self.status_code = 200
        self.content = text.encode("utf-8")


class _FakeSession:
    def __init__(self, text):
        self._text = text

    def get(self, *a, **k):
        return _FakeResp(self._text)


old_session = tr.tx.HTTP_SESSION
tr.tx.HTTP_SESSION = _FakeSession(trrust_text)
try:
    fetched = tr.fetch_trrust_regulons("9606", force_refresh=True)
finally:
    tr.tx.HTTP_SESSION = old_session

tf1 = fetched["regulons"].get("TF1")
check(tf1 is not None, "TF1 parsed from TRRUST TSV")
if tf1:
    check(tf1["genes"] == {"G1", "G2", "G3"}, f"TF1 genes includes ALL targets incl. Unknown-mode (got {tf1['genes']})")
    check(tf1["weights"] == {"G1": 1, "G2": -1}, f"TF1 weights EXCLUDES the Unknown-mode target (got {tf1['weights']})")
tf2 = fetched["regulons"].get("TF2")
check(tf2 is not None and tf2["weights"] == {"G1": 1}, "TF2 parsed correctly (single Activation edge)")


print("=== ChEA3 GMT parsing (cell-line union + no-description-column format) ===")
chea3_gmt_text = (
    "STAT1_HELAS3_HG19\tG1\tG2\tG3\n"
    "STAT1_K562_HG19\tG3\tG4\n"
    "JUN_GM12878_HG19\tG5\tG6\n"
)


class _FakeChea3Resp:
    def __init__(self, text):
        self.status_code = 200
        self.content = text.encode("utf-8")


class _FakeChea3Session:
    def get(self, *a, **k):
        return _FakeChea3Resp(chea3_gmt_text)


old_session = tr.tx.HTTP_SESSION
tr.tx.HTTP_SESSION = _FakeChea3Session()
try:
    fetched_chea3 = tr.fetch_chea3_regulons(libraries=["ENCODE_ChIP-seq"], force_refresh=True)
finally:
    tr.tx.HTTP_SESSION = old_session

stat1 = fetched_chea3["regulons"].get("STAT1")
check(stat1 is not None, "STAT1 parsed and merged across its two cell-line rows")
if stat1:
    check(stat1["genes"] == {"G1", "G2", "G3", "G4"}, f"STAT1 targets UNIONED across cell lines (got {stat1['genes']})")
    check(stat1["weights"] == {}, "STAT1 weights empty (ChEA3 is unsigned)")
check(fetched_chea3["signed"] is False, "ChEA3 fetch result flagged signed=False")
jun = fetched_chea3["regulons"].get("JUN")
check(jun is not None and jun["genes"] == {"G5", "G6"}, "JUN parsed from its own single row")


print("=== OmniPath JSON row parsing (sign fields + COMPLEX: exclusion) ===")
omnipath_rows = [
    {"source": "P01", "target": "P02", "source_genesymbol": "TF3", "target_genesymbol": "G1",
     "is_stimulation": True, "is_inhibition": False},
    {"source": "P01", "target": "P03", "source_genesymbol": "TF3", "target_genesymbol": "G2",
     "is_stimulation": False, "is_inhibition": True},
    {"source": "P01", "target": "P04", "source_genesymbol": "TF3", "target_genesymbol": "G3",
     "is_stimulation": True, "is_inhibition": True},  # ambiguous -> genes only, no weight
    {"source": "COMPLEX:P01_P05", "target": "P06", "source_genesymbol": "TF3_TF4", "target_genesymbol": "G4",
     "is_stimulation": True, "is_inhibition": False},  # must be SKIPPED entirely
]
built = tr._build_regulons_from_omnipath_rows(omnipath_rows)
check("TF3_TF4" not in built, "protein-complex source row excluded entirely (COMPLEX: prefix)")
tf3 = built.get("TF3")
check(tf3 is not None, "TF3 built from OmniPath rows")
if tf3:
    check(tf3["genes"] == {"G1", "G2", "G3"}, f"TF3 genes includes the ambiguous edge too (got {tf3['genes']})")
    check(tf3["weights"] == {"G1": 1, "G2": -1}, f"TF3 weights excludes the ambiguous (both-true) edge (got {tf3['weights']})")


# =============================================================================
# 6. Cache round-trip
# =============================================================================
print("=== cache serialize/deserialize round-trip ===")
original_regulons = {
    "TFX": {"name": "TFX", "genes": {"G1", "G2", "G3"}, "weights": {"G1": 1, "G2": -1}},
}
serialized = tr._serialize_regulons(original_regulons)
check(isinstance(serialized["TFX"]["genes"], list), "genes serialized as a list (JSON-safe)")
deserialized = tr._deserialize_regulons(serialized)
check(deserialized["TFX"]["genes"] == original_regulons["TFX"]["genes"], "genes round-trip as an equal set")
check(deserialized["TFX"]["weights"] == original_regulons["TFX"]["weights"], "weights round-trip exactly")


# =============================================================================
# 7. run_regulon_enrichment: orchestration around the reused _run_ora_for_sets
#    (the hypergeometric math itself is already covered by the rest of the
#    suite via Genomics/Transcriptomics ORA - this checks THIS module's own
#    wiring: database/direction looping, fetch-failure messaging, and that an
#    UNSIGNED database like ChEA3 participates normally here, unlike activity).
# =============================================================================
print("=== run_regulon_enrichment: orchestration ===")

enrichment_resolved = _make_resolved_df({
    **{f"DEG{i}": 3.0 for i in range(6)},        # up-regulated, in TF_ENRICH's regulon
    **{f"BG{i}": 0.01 for i in range(24)},        # background, not in any regulon
}, pvalue=0.001)

enrich_fake = {
    "regulons": {
        "TF_ENRICH": {"name": "TF_ENRICH", "genes": {f"DEG{i}" for i in range(6)}, "weights": {}},
        "TF_OTHER": {"name": "TF_OTHER", "genes": {f"BG{i}" for i in range(5)}, "weights": {}},
    },
    "signed": False, "fetched_at": "x", "fetch_failed": False, "note": None,
}
original_fetch = _patch_fetch(enrich_fake)
try:
    enrich_result = tr.run_regulon_enrichment(
        enrichment_resolved, ["ChEA3"], "9606",
        padj_cutoff=0.05, lfc_cutoff=1.0, min_term_size=3, max_term_size=None,
    )
finally:
    tr.fetch_regulons = original_fetch

up_table = enrich_result["tables"]["ChEA3 (up)"]["table"]
check(len(up_table) >= 1, f"ChEA3 (up) produced at least one tested term (got {len(up_table)} rows)")
if len(up_table) >= 1:
    row = up_table[up_table["Term ID"] == "TF_ENRICH"].iloc[0]
    check(int(row["k"]) == 6, f"TF_ENRICH overlap k == 6 (all 6 DEGs are its targets) (got {row['k']})")
    check(row["p-value"] < 0.01, f"TF_ENRICH is strongly enriched (got p={row['p-value']})")

# Fetch failure must produce an explanatory empty block, not raise or silently
# vanish the database from the results.
failed_fetch = {"regulons": {}, "signed": True, "fetched_at": "x", "fetch_failed": True, "note": "boom"}
original_fetch = _patch_fetch(failed_fetch)
try:
    failed_result = tr.run_regulon_enrichment(
        enrichment_resolved, ["TRRUST"], "9606", min_term_size=3, max_term_size=None,
    )
finally:
    tr.fetch_regulons = original_fetch
failed_block = failed_result["tables"]["TRRUST (up)"]
check(failed_block["table"].empty, "fetch-failed database yields an empty (not missing) table")
check("Could not fetch" in failed_block["note"], f"fetch-failed note explains why (got: {failed_block['note']})")


# =============================================================================
# 7b. run_regulon_activity: tables dict preserves the CALLER's database order,
#     not "unsigned first" - the client's results screen shows
#     Object.keys(tables)[0] by default, so an unsigned database silently
#     winning that slot means the default view is always empty regardless of
#     what the user actually selected. Caught live via Playwright: with
#     databases=["TRRUST","ChEA3"], the UI opened on ChEA3 (0 TFs) every time.
# =============================================================================
print("=== run_regulon_activity: tables dict preserves caller order (ChEA3-first regression) ===")
original_fetch = _patch_fetch(fake_regulons)  # TF1 fixture from earlier, signed
try:
    ordered_result = tr.run_regulon_activity(resolved_up, ["TRRUST", "ChEA3"], "9606", metric="log2fc", min_regulon_size=3)
finally:
    tr.fetch_regulons = original_fetch
ordered_keys = list(ordered_result["tables"].keys())
check(ordered_keys == ["TRRUST", "ChEA3"],
      f"tables dict key order matches the caller's database order, not unsigned-first (got {ordered_keys})")

reversed_result_fetch = _patch_fetch(fake_regulons)
try:
    reversed_result = tr.run_regulon_activity(resolved_up, ["ChEA3", "TRRUST"], "9606", metric="log2fc", min_regulon_size=3)
finally:
    tr.fetch_regulons = reversed_result_fetch
reversed_keys = list(reversed_result["tables"].keys())
check(reversed_keys == ["ChEA3", "TRRUST"], f"order tracks the caller's list even when unsigned is listed first (got {reversed_keys})")


# =============================================================================
# 8. ChEA3 human-only gate (fetch_regulons dispatcher, not fetch_chea3_regulons
#    itself - mirrors run_ora's MSigDB Hallmarks human-only gate)
# =============================================================================
print("=== ChEA3 human-only gate ===")
mouse_chea3 = tr.fetch_regulons("ChEA3", "10090")
check(not mouse_chea3["regulons"], "ChEA3 returns no regulons for a non-human taxon")
check(not mouse_chea3["fetch_failed"], "gated (not failed) - this is an expected, explained skip")
check(mouse_chea3.get("note") and "human" in mouse_chea3["note"].lower(),
      f"note explains the human-only restriction (got: {mouse_chea3.get('note')})")


# =============================================================================
# 9. Permutation null (_null_slopes_shared / _permutation_p_value) - this is
#    the actual point of this session's work: replacing an anti-conservative
#    parametric p-value with one that makes no independence assumption.
# =============================================================================
print("=== _null_slopes_shared: null-sharing correctness ===")
rng_a = np.random.default_rng(123)
big_n = 500
big_y = np.random.default_rng(7).normal(size=big_n)
big_yc = big_y - big_y.mean()
null_map = tr._null_slopes_shared([(5, 3), (5, 3), (10, 2)], big_yc, big_n, 200, rng_a)
check(set(null_map.keys()) == {(5, 3), (10, 2)}, f"duplicate signature deduplicated (got {set(null_map.keys())})")
check(null_map[(5, 3)].shape == (200,), f"null array has n_permutations entries (got {null_map[(5, 3)].shape})")
check(not np.array_equal(null_map[(5, 3)], null_map[(10, 2)]),
      "two DIFFERENT signatures get different null arrays (not accidentally aliased)")

print("=== _null_slopes_shared: determinism (same seed -> identical result) ===")
rng_b1 = np.random.default_rng(99)
rng_b2 = np.random.default_rng(99)
null_1 = tr._null_slopes_shared([(4, 4)], big_yc, big_n, 150, rng_b1)
null_2 = tr._null_slopes_shared([(4, 4)], big_yc, big_n, 150, rng_b2)
check(np.array_equal(null_1[(4, 4)], null_2[(4, 4)]), "identical seed produces identical null array")

print("=== monotonicity: ranking by |beta1| matches ranking by |t| at a fixed signature ===")
# Build several regulons of the SAME signature (5,3) on the same big_y, drawn
# from random disjoint index sets, and confirm sorting by |beta1| and by |t|
# (the parametric score) gives the identical order - this is the assumption
# _null_slopes_shared's docstring rests the "store slopes only" choice on.
rng_mono = np.random.default_rng(55)
mono_rows = []
big_syy = float(np.dot(big_yc, big_yc))
for _ in range(15):
    idx = rng_mono.choice(big_n, size=8, replace=False)
    pos_idx, neg_idx = list(idx[:5]), list(idx[5:])
    sc, pv, nu, b1 = tr._ulm_score(pos_idx, neg_idx, big_yc, big_syy, big_n)
    if sc is not None:
        mono_rows.append((abs(sc), abs(b1)))
order_by_score = sorted(range(len(mono_rows)), key=lambda i: mono_rows[i][0])
order_by_beta1 = sorted(range(len(mono_rows)), key=lambda i: mono_rows[i][1])
check(order_by_score == order_by_beta1, "ranking by |t| and by |beta1| agree at a fixed (n_pos, n_neg) signature")

print("=== permutation vs parametric: agree under INDEPENDENT (iid) gene values ===")
rng_iid = np.random.default_rng(2024)
n_iid = 400
y_iid = rng_iid.normal(size=n_iid)
yc_iid = y_iid - y_iid.mean()
syy_iid = float(np.dot(yc_iid, yc_iid))
# A TF regulon with genuinely random (not co-regulated) targets - the t-test's
# independence assumption HOLDS here, so permutation and parametric p should
# be close.
idx_iid = rng_iid.choice(n_iid, size=20, replace=False)
pos_iid, neg_iid = list(idx_iid[:12]), list(idx_iid[12:])
sc_iid, p_param_iid, _, b1_iid = tr._ulm_score(pos_iid, neg_iid, yc_iid, syy_iid, n_iid)
null_iid = tr._null_slopes_shared([(12, 8)], yc_iid, n_iid, 3000, np.random.default_rng(11))[(12, 8)]
p_perm_iid = tr._permutation_p_value(abs(b1_iid), null_iid)
check(abs(p_perm_iid - p_param_iid) < 0.03,
      f"under iid gene values, permutation p and parametric p agree closely (perm={p_perm_iid:.4f}, param={p_param_iid:.4f})")

print("=== permutation vs parametric under correlated (co-regulated) targets: NEITHER is asserted superior ===")
# IMPORTANT FINDING (checked directly, not assumed): an earlier version of
# this test asserted permutation must come out MATERIALLY MORE CONSERVATIVE
# than the parametric p-value under target correlation - i.e. that
# permutation "fixes" the parametric test's independence-assumption weakness.
# That assertion was FALSE and is deliberately not made here. Gene-label
# permutation draws its null from RANDOM same-sized subsets of the ranked
# list, which are typically uncorrelated - so it can UNDER-represent how
# extreme a genuinely co-regulated regulon's statistic is, in the SAME
# direction as the parametric test's own weakness. This mirrors the already-
# documented limitation of gene-label permutation for GSEA itself (see
# run_gsea's module docstring in Transcriptomics_module.py). This test only
# checks that BOTH p-values are computed validly (in range, and genuinely
# different NUMBERS from two different mechanisms) - it does not, and should
# not, assert a direction of disagreement.
rng_corr = np.random.default_rng(4242)
n_corr = 400
block_shift = rng_corr.normal(scale=2.5)  # one shared shift for the whole block
y_corr = rng_corr.normal(size=n_corr)
y_corr[:50] += block_shift + rng_corr.normal(scale=0.3, size=50)  # correlated block, small within-block noise
yc_corr = y_corr - y_corr.mean()
syy_corr = float(np.dot(yc_corr, yc_corr))
pos_corr = list(range(12))       # entirely inside the correlated block
neg_corr = list(range(12, 20))   # outside it - independent noise
sc_corr, p_param_corr, _, b1_corr = tr._ulm_score(pos_corr, neg_corr, yc_corr, syy_corr, n_corr)
null_corr = tr._null_slopes_shared([(12, 8)], yc_corr, n_corr, 3000, np.random.default_rng(13))[(12, 8)]
p_perm_corr = tr._permutation_p_value(abs(b1_corr), null_corr)
check(0.0 < p_param_corr <= 1.0, f"parametric p-value valid under correlation (got {p_param_corr})")
check(0.0 < p_perm_corr <= 1.0, f"permutation p-value valid under correlation (got {p_perm_corr})")

print("=== _thin_null round-trip: thinned quantiles track the full null ===")
full_null = tr._null_slopes_shared([(6, 4)], big_yc, big_n, 2000, np.random.default_rng(31))[(6, 4)]
thinned = tr.tx._thin_null(full_null, target=200)
check(len(thinned) <= 200, f"thinned null respects the target cap (got {len(thinned)})")
check(abs(np.median(thinned) - np.median(full_null)) < 0.05 * (np.std(full_null) + 1e-9),
      "thinned null's median tracks the full null's median")
check(abs(min(thinned) - float(full_null.min())) < 1e-3, "thinned null preserves the minimum (order statistics, within the function's own 4dp rounding)")
check(abs(max(thinned) - float(full_null.max())) < 1e-3, "thinned null preserves the maximum (order statistics, within the function's own 4dp rounding)")

print("=== run_regulon_activity: end-to-end permutation columns present and internally consistent ===")
result_perm = None
original_fetch = _patch_fetch(fake_regulons)
try:
    result_perm = tr.run_regulon_activity(
        resolved_up, ["TRRUST"], "9606", metric="log2fc", min_regulon_size=3, n_permutations=500, seed=7,
    )
finally:
    tr.fetch_regulons = original_fetch
tbl = result_perm["tables"]["TRRUST"]["table"]
check(list(tbl.columns) == tr.ACTIVITY_RESULT_COLUMNS, f"table columns match ACTIVITY_RESULT_COLUMNS (got {list(tbl.columns)})")
check("p-value (permutation)" in tbl.columns and "FDR (permutation)" in tbl.columns, "permutation columns present")
row0 = tbl.iloc[0]
check(0.0 < row0["p-value (permutation)"] <= 1.0, f"permutation p-value in (0,1] (got {row0['p-value (permutation)']})")
detail = result_perm["tables"]["TRRUST"]["tf_regression_detail"]
check("TF1" in detail, "top TF has a stored regression detail entry")
if "TF1" in detail:
    sig_key = detail["TF1"]["signature"]
    check(sig_key in result_perm["tables"]["TRRUST"]["null_slopes_by_signature"],
          "TF1's signature key resolves into the stored null_slopes_by_signature dict")
    check(set(detail["TF1"]["activated_symbols"]) == {"A1", "A2", "A3"}, "activated target symbols stored correctly")
    check(set(detail["TF1"]["repressed_symbols"]) == {"R1", "R2"}, "repressed target symbols stored correctly")

# Determinism at the run_regulon_activity level: same seed -> identical
# permutation p-values across two independent calls.
original_fetch = _patch_fetch(fake_regulons)
try:
    result_perm_2 = tr.run_regulon_activity(
        resolved_up, ["TRRUST"], "9606", metric="log2fc", min_regulon_size=3, n_permutations=500, seed=7,
    )
finally:
    tr.fetch_regulons = original_fetch
p1 = result_perm["tables"]["TRRUST"]["table"].iloc[0]["p-value (permutation)"]
p2 = result_perm_2["tables"]["TRRUST"]["table"].iloc[0]["p-value (permutation)"]
check(p1 == p2, f"same seed reproduces the identical permutation p-value across runs (got {p1} vs {p2})")


# =============================================================================
# 10. Sensitivity check: deliberately break _null_slopes_shared's signature
#     handling and confirm the null-sharing test goes red (established rigor
#     pattern this session - a test that cannot fail is not evidence).
# =============================================================================
print("=== sensitivity check: null-sharing test is actually discriminating ===")
_original_null_slopes_shared = tr._null_slopes_shared


def _broken_null_slopes_shared(signatures, yc, n, n_permutations, rng, chunk=200):
    # DELIBERATELY BROKEN: collapses every signature onto the SAME array,
    # simulating the bug this test exists to catch (a signature-keying error
    # that aliases two different TFs' nulls together).
    real = _original_null_slopes_shared(signatures, yc, n, n_permutations, rng, chunk)
    keys = list(real.keys())
    if len(keys) >= 2:
        for k in keys[1:]:
            real[k] = real[keys[0]]
    return real


tr._null_slopes_shared = _broken_null_slopes_shared
broken_null_map = tr._null_slopes_shared([(5, 3), (10, 2)], big_yc, big_n, 200, np.random.default_rng(123))
sensitivity_caught = np.array_equal(broken_null_map[(5, 3)], broken_null_map[(10, 2)])
tr._null_slopes_shared = _original_null_slopes_shared
check(sensitivity_caught, "deliberately-aliased signatures WOULD have been caught (the real test asserts inequality)")


# =============================================================================
# 11. Regression-figure payload correctness. The Regulon Detail figure claims
#     to draw "the actual regression". An earlier version did NOT: it plotted
#     the client's log2FC lookup while the fit ran on the RANKING METRIC
#     (default signed -log10(p)), and pivoted its lines through a centroid
#     computed over the TARGETS instead of over all tested genes - measured
#     141x off on real data. The Playwright test at the time only asserted
#     "a trace exists", so it passed anyway. These checks assert the payload
#     is self-consistent enough for the figure to be drawable correctly.
# =============================================================================
print("=== regression-figure payload: values are the RANKING METRIC, not log2FC ===")
# Fixture where log2FC and signed_log_pvalue are deliberately DIFFERENT, so a
# regression to the old (log2FC) behaviour is detectable rather than masked.
distinct_genes = {"A1": 3.0, "A2": 2.5, "A3": 2.8, "R1": -2.5, "R2": -3.0}
for i in range(15):
    distinct_genes[f"FILLER{i}"] = 0.05 * (-1 if i % 2 else 1)
resolved_distinct = _make_resolved_df(distinct_genes, pvalue=0.001)   # p far from 10^-|log2FC|

original_fetch = _patch_fetch(fake_regulons)
try:
    fig_result = tr.run_regulon_activity(
        resolved_distinct, ["TRRUST"], "9606",
        metric="signed_log_pvalue", min_regulon_size=3, n_permutations=200, seed=5,
    )
finally:
    tr.fetch_regulons = original_fetch

fig_block = fig_result["tables"]["TRRUST"]
fig_detail = fig_block["tf_regression_detail"]["TF1"]

check("activated_values" in fig_detail and "repressed_values" in fig_detail,
      "detail ships per-target y-VALUES, not just symbols")
check(len(fig_detail["activated_values"]) == len(fig_detail["activated_symbols"]),
      "activated_values is aligned 1:1 with activated_symbols")
check(len(fig_detail["repressed_values"]) == len(fig_detail["repressed_symbols"]),
      "repressed_values is aligned 1:1 with repressed_symbols")

# The shipped values must equal build_ranked_list's metric values, NOT log2FC.
ranked_check = tr.tx.build_ranked_list(resolved_distinct, "signed_log_pvalue")
metric_lookup = dict(zip(ranked_check["symbols"], ranked_check["values"]))
log2fc_lookup = {str(r["gene_symbol"]).upper(): float(r["log2FC"]) for _, r in resolved_distinct.iterrows()}
shipped = dict(zip(fig_detail["activated_symbols"], fig_detail["activated_values"]))
matches_metric = all(abs(shipped[s] - metric_lookup[s]) < 1e-3 for s in shipped)
matches_log2fc = all(abs(shipped[s] - log2fc_lookup[s]) < 1e-3 for s in shipped)
check(matches_metric, f"shipped y-values match the RANKING METRIC the fit used (got {shipped})")
check(not matches_log2fc, "shipped y-values are NOT log2FC (the old, wrong source) - fixture makes these distinguishable")

print("=== regression-figure payload: centroid matches the fit's own, not the targets' ===")
n_genes_fig = fig_block["n_genes"]
n_pos_fig = len(fig_detail["activated_symbols"])
n_neg_fig = len(fig_detail["repressed_symbols"])
expected_x_mean = (n_pos_fig - n_neg_fig) / n_genes_fig
check(abs(fig_detail["x_mean"] - expected_x_mean) < 1e-12,
      f"x_mean == (n_pos-n_neg)/n over ALL tested genes (got {fig_detail['x_mean']}, expected {expected_x_mean})")
wrong_x_mean = (n_pos_fig - n_neg_fig) / (n_pos_fig + n_neg_fig)
check(abs(fig_detail["x_mean"] - wrong_x_mean) > 1e-6,
      "x_mean is NOT the over-the-targets value the old figure used (they are distinguishable here)")
check(abs(fig_block["y_mean"] - float(np.mean(ranked_check["values"]))) < 1e-9,
      "y_mean is the mean ranking value over ALL tested genes")

print("=== regression-figure payload: the drawn line IS the OLS fit (brute-force check) ===")
# NOTE ON WHAT IS *NOT* ASSERTED HERE. An earlier version of this test
# demanded the fitted line land inside the target point cloud. That was
# measured against real data and found to be FALSE for legitimate reasons:
# with x = {+1, -1, 0} and the x=0 group holding ~99% of genes, the intercept
# is pinned by that mass, so a lopsided regulon (real example: RBMX, 1
# activated vs 6 repressed) correctly produces a line that passes nowhere
# near its single activated point. "Lands inside the points" was never a
# property of this model - it was my assumption. The correct check is that
# the shipped (beta1, x_mean, y_mean) reproduce a brute-force full-length
# OLS, which is what the figure draws from.
xs = np.zeros(fig_block["n_genes"])
sym_pos = {s: i for i, s in enumerate(ranked_check["symbols"])}
for s in fig_detail["activated_symbols"]:
    xs[sym_pos[s]] = 1.0
for s in fig_detail["repressed_symbols"]:
    xs[sym_pos[s]] = -1.0
ys = np.asarray(ranked_check["values"], dtype=float)
xm_bf, ym_bf = xs.mean(), ys.mean()
b1_bf = float(np.dot(xs - xm_bf, ys - ym_bf) / np.dot(xs - xm_bf, xs - xm_bf))
check(abs(fig_detail["beta1"] - b1_bf) < 1e-9,
      f"shipped beta1 reproduces brute-force full-length OLS (got {fig_detail['beta1']:.6f}, ref {b1_bf:.6f})")
check(abs(fig_detail["x_mean"] - xm_bf) < 1e-12, "shipped x_mean reproduces brute-force OLS x-centroid")
check(abs(fig_block["y_mean"] - ym_bf) < 1e-9, "shipped y_mean reproduces brute-force OLS y-centroid")
# The centroid form the figure uses must be algebraically identical to the
# intercept form - i.e. drawing y = y_mean + beta1*(x - x_mean) really is
# drawing y = beta0 + beta1*x.
beta0_bf = ym_bf - b1_bf * xm_bf
centroid_form = fig_block["y_mean"] + fig_detail["beta1"] * (1 - fig_detail["x_mean"])
intercept_form = beta0_bf + b1_bf * 1.0
check(abs(centroid_form - intercept_form) < 1e-9,
      "the figure's centroid line form is identical to beta0 + beta1*x")

print("=== regression-figure payload: all-genes distribution sample for the x=0 violin ===")
check(len(fig_block["ranking_values_sample"]) > 0, "ranking_values_sample is populated")
check(len(fig_block["ranking_values_sample"]) <= 400, "ranking_values_sample respects its cap")
check(fig_block["metric"] == "signed_log_pvalue", "block carries the metric name for the axis label")


# =============================================================================
# 12. Empirical-null calibration (_empirical_null_calibration). Built ONLY
#     after tests/calibration_regulon_activity.py measured, on real
#     data, a material inflation (lambda ~1.45-1.67, 22-36% false-positive
#     rate at nominal alpha=0.05) - see that script and this function's own
#     docstring for the full justification, including the exchangeability
#     premise these tests check rather than assume.
# =============================================================================
print("=== _empirical_null_calibration: no-op when there is nothing to calibrate ===")
rng_noop = np.random.default_rng(11)
# Independent (iid) scores - by construction, no inflation should be detected.
iid_scores = rng_noop.normal(loc=0.0, scale=1.0, size=200)
center_noop, lam_noop, n_used_noop, note_noop = tr._empirical_null_calibration(iid_scores)
check(n_used_noop == 200, f"n_used == 200 (got {n_used_noop})")
check(note_noop is None, f"no skip-note when calibration actually runs (got {note_noop!r})")
check(abs(lam_noop - 1.0) < 0.15, f"lambda stays close to 1.0 on genuinely independent scores (got {lam_noop:.4f})")

print("=== _empirical_null_calibration: bites when target co-regulation inflates the spread ===")
# Scores drawn from a WIDER distribution (sigma=1.8) simulate the effect
# co-regulation has on the observed score spread - the estimator must
# recover that inflation ON AVERAGE. Averaged over reps rather than a single
# draw: the truncated-normal estimator (see section 12b below) restricts
# itself to a central window, so at n=200 a SINGLE draw has real sampling
# variance (empirically std ~0.2 at true sigma=1.8) even though it is
# unbiased across repetitions (measured mean 1.76 at true 1.8, n=200) - a
# single-draw threshold tight enough to be meaningful would be flaky.
inflated_runs = []
inflated_notes = []
for rep in range(30):
    rng_inf_rep = np.random.default_rng(3000 + rep)
    inflated_runs.append(tr._empirical_null_calibration(rng_inf_rep.normal(loc=0.0, scale=1.8, size=200)))
lam_inf_mean = float(np.mean([r[1] for r in inflated_runs]))
check(lam_inf_mean > 1.5, f"lambda recovers the true inflation factor on average (~1.8 expected, got mean {lam_inf_mean:.4f} over 30 reps)")
check(all(r[3] is None for r in inflated_runs), "no skip-note when calibration genuinely applies")

print("=== _empirical_null_calibration: the clamp holds - never below 1.0 ===")
# A distribution NARROWER than theory (sigma=0.5): dividing by a sub-1.0
# sigma_hat would make results LESS conservative, which the clamp forbids.
narrow_scores = rng_noop.normal(loc=0.0, scale=0.5, size=200)
center_narrow, lam_narrow, n_used_narrow, note_narrow = tr._empirical_null_calibration(narrow_scores)
check(lam_narrow == 1.0, f"lambda clamped to exactly 1.0 for a narrower-than-theory distribution (got {lam_narrow})")

print("=== _empirical_null_calibration: guards ===")
_, lam_few, n_few, note_few = tr._empirical_null_calibration([0.1, 0.5, -0.3])
check(lam_few == 1.0, "too few TFs -> lambda left at 1.0")
check(note_few is not None and "3" in note_few, f"too-few-TFs note explains why (got: {note_few!r})")

_, lam_degenerate, n_deg, note_degenerate = tr._empirical_null_calibration([1.0] * 50)
check(lam_degenerate == 1.0, "zero-MAD (all-identical) scores -> lambda left at 1.0")
check(note_degenerate is not None and "zero spread" in note_degenerate.lower(),
      f"zero-MAD note explains why (got: {note_degenerate!r})")

print("=== _empirical_null_calibration: sensitivity - deliberately break the clamp, confirm red ===")
_original_calibration = tr._empirical_null_calibration


def _broken_calibration_no_clamp(scores, min_tfs=tr.REGULON_MIN_TFS_FOR_CALIBRATION):
    # DELIBERATELY BROKEN: removes the >= 1.0 clamp entirely.
    arr = np.asarray(scores, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size < min_tfs:
        return 0.0, 1.0, int(arr.size), "skipped"
    center = float(np.median(arr))
    mad = float(np.median(np.abs(arr - center)))
    if mad <= 0:
        return center, 1.0, int(arr.size), "skipped"
    return center, mad / 0.6745, int(arr.size), None


tr._empirical_null_calibration = _broken_calibration_no_clamp
_, lam_broken, _, _ = tr._empirical_null_calibration(narrow_scores)
clamp_test_would_catch_it = lam_broken != 1.0
tr._empirical_null_calibration = _original_calibration
check(clamp_test_would_catch_it, "the clamp test IS sensitive: an unclamped implementation gives lambda < 1.0 here, which the real test rejects")


# =============================================================================
# 12b. Truncated-normal null-scale estimator (_truncated_normal_sigma),
#      promoted from tests/probe_better_null_estimator.py and
#      probe_truncated_validity.py after those probes justified REPLACING
#      plain MAD/0.6745 as the estimator _empirical_null_calibration uses.
#      _empirical_null_calibration itself now RETURNS the truncated lambda -
#      "plain MAD" below is a deliberately separate, local reimplementation of
#      the OLD pre-fix formula, kept only as the baseline these tests measure
#      the improvement against.
# =============================================================================
def _plain_mad_lambda(scores):
    """The estimator this module used BEFORE the truncated-normal MLE
    replaced it - reimplemented locally (not imported) because the module no
    longer exposes it standalone; this is the baseline, not the code under
    test."""
    arr = np.asarray(scores, dtype=float)
    arr = arr[np.isfinite(arr)]
    center = float(np.median(arr))
    mad = float(np.median(np.abs(arr - center)))
    if mad <= 0:
        return 1.0
    return max(1.0, mad / 0.6745)


print("=== truncated estimator vs plain MAD: fooled-by-signal grid (truth = 1.0 in every row) ===")
# At the TF counts this module actually produces (TRRUST ~134, CollecTRI
# ~494). Any lambda > 1.0 here is the estimator being fooled by real signal,
# not detecting genuine inflation - there is none in this fixture.
for n_tfs in (134, 494):
    for frac in (0.10, 0.25, 0.40):
        mad_runs, trunc_runs = [], []
        for rep in range(40):
            rng_grid = np.random.default_rng(1000 + rep)
            s = rng_grid.normal(size=n_tfs)
            n_act = int(frac * n_tfs)
            if n_act:
                s[:n_act] += rng_grid.choice([-1, 1], size=n_act) * rng_grid.uniform(4, 9, size=n_act)
            _, lam_trunc, _, _ = tr._empirical_null_calibration(s, min_tfs=10)
            lam_mad = _plain_mad_lambda(s)
            trunc_runs.append(lam_trunc)
            mad_runs.append(lam_mad)
        mad_m = float(np.mean(mad_runs))
        trunc_m = float(np.mean(trunc_runs))
        check(trunc_m < mad_m - 0.03,
              f"n_tfs={n_tfs}, {int(frac*100)}% active: truncated lambda ({trunc_m:.3f}) beats plain MAD "
              f"({mad_m:.3f}) by more than 0.03")

print("=== truncated estimator: still recovers GENUINE inflation (the condition that matters most) ===")
# A "fix" that flattened lambda to 1.0 would silently destroy calibration -
# this is the case that would catch that regression. No signal at all here;
# the whole null is genuinely scaled by true_sigma, and the estimator must
# report close to that, not close to 1.0.
for true_sigma in (1.3, 1.6, 2.0, 3.0):
    trunc_runs = []
    for rep in range(40):
        rng_inf = np.random.default_rng(500 + rep)
        s = rng_inf.normal(scale=true_sigma, size=300)
        _, lam_trunc, _, _ = tr._empirical_null_calibration(s, min_tfs=10)
        trunc_runs.append(lam_trunc)
    trunc_m = float(np.mean(trunc_runs))
    rel_err = abs(trunc_m - true_sigma) / true_sigma
    check(rel_err < 0.08,
          f"true sigma={true_sigma}: truncated lambda ({trunc_m:.3f}) recovers it within 8% (rel err {rel_err:.1%})")

print("=== truncated estimator: graded (continuous) signal - materially better, not perfect (kept separate) ===")
# Weaker claim than the clump case above, deliberately not folded into it:
# with continuous effect sizes some genuinely-active TFs land INSIDE the
# truncation window, so the correction is only partial. Asserts "better than
# plain MAD", not "close to 1.0" - a stronger assertion here would be false
# and would have to be walked back, exactly the kind of overclaim this
# module's docstrings were rewritten to avoid.
for frac in (0.10, 0.25, 0.40):
    mad_runs, trunc_runs = [], []
    for rep in range(40):
        rng_graded = np.random.default_rng(900 + rep)
        s = rng_graded.normal(size=300)
        n_act = int(frac * 300)
        if n_act:
            s[:n_act] += rng_graded.normal(scale=3.0, size=n_act)
        _, lam_trunc, _, _ = tr._empirical_null_calibration(s, min_tfs=10)
        lam_mad = _plain_mad_lambda(s)
        trunc_runs.append(lam_trunc)
        mad_runs.append(lam_mad)
    mad_m = float(np.mean(mad_runs))
    trunc_m = float(np.mean(trunc_runs))
    check(trunc_m < mad_m - 0.03,
          f"graded signal, {int(frac*100)}% active: truncated lambda ({trunc_m:.3f}) still beats plain MAD "
          f"({mad_m:.3f}), even though neither is close to the 1.0 truth here")

print("=== truncated estimator: no calibrated p-value is ever exactly 0 after the floor ===")
extreme_scores = np.concatenate([np.array([50.0, -50.0]), rng_noop.normal(size=150)])
_, lam_extreme, _, _ = tr._empirical_null_calibration(extreme_scores, min_tfs=10)
from scipy.special import ndtr as _ndtr_check
z_extreme = 50.0 / lam_extreme
p_extreme = max(2.0 * float(_ndtr_check(-abs(z_extreme))), tr.tx.GSEA_MIN_P)
check(p_extreme > 0.0, f"an extreme |z| still yields a strictly positive calibrated p-value (got {p_extreme:.3e})")
check(p_extreme >= tr.tx.GSEA_MIN_P, f"and it is never below the floor tr.tx.GSEA_MIN_P (got {p_extreme:.3e} vs floor {tr.tx.GSEA_MIN_P:.3e})")

print("=== sensitivity: deliberately break the truncation bias correction, confirm inflation-recovery goes red ===")
# Drops the "factor" term entirely (as if the central window needed no
# correction for having thrown the tails away) - i.e. plain central-window
# variance. This should UNDER-estimate true inflation, since a truncated
# window's raw variance is smaller than the untruncated population variance
# it's supposed to stand in for.
_source_before = None
with open(tr.__file__, "rb") as _f:
    _source_before = _f.read()
_original_truncated_sigma = tr._truncated_normal_sigma


def _broken_truncated_sigma(deviations, sigma0, c_mult=1.5, iters=60):
    sigma = sigma0
    for _ in range(iters):
        c = c_mult * sigma
        inner = deviations[np.abs(deviations) < c]
        if inner.size < 10:
            return sigma0
        m2 = float(np.mean(inner ** 2))
        # BROKEN: no truncation-bias correction (factor forced to 1.0).
        new_sigma = math.sqrt(m2)
        if not math.isfinite(new_sigma) or new_sigma <= 0:
            return sigma0
        if abs(new_sigma - sigma) < 1e-9:
            return new_sigma
        sigma = new_sigma
    return sigma


tr._truncated_normal_sigma = _broken_truncated_sigma
broken_rel_errs = []
for true_sigma in (1.3, 1.6, 2.0, 3.0):
    trunc_runs = []
    for rep in range(40):
        rng_broken = np.random.default_rng(500 + rep)
        s = rng_broken.normal(scale=true_sigma, size=300)
        _, lam_broken_trunc, _, _ = tr._empirical_null_calibration(s, min_tfs=10)
        trunc_runs.append(lam_broken_trunc)
    trunc_m = float(np.mean(trunc_runs))
    broken_rel_errs.append(abs(trunc_m - true_sigma) / true_sigma)
tr._truncated_normal_sigma = _original_truncated_sigma
with open(tr.__file__, "rb") as _f:
    _source_after = _f.read()
check(_source_after == _source_before,
      "the module's source file on disk was never touched by this sensitivity test (monkeypatch only, in-memory)")
check(any(e >= 0.08 for e in broken_rel_errs),
      f"with the bias correction removed, at least one true-sigma case now fails the 8% recovery tolerance that "
      f"passed above (rel errs: {[f'{e:.1%}' for e in broken_rel_errs]}) - proves the test above is not vacuous")


print("=== exchangeability premise: the t-statistic's null does not depend on regulon size/split ===")
# The whole pooling argument rests on every TF's t-statistic sharing ONE null
# distribution under H0, regardless of (n_pos, n_neg). Verify empirically:
# simulate many null (no-true-effect) TFs across DELIBERATELY varied
# signatures and confirm the resulting t-statistics have equal variance
# across strata (an F-test-free, distribution-free check: compare stratum
# standard deviations directly against the pooled one).
rng_exch = np.random.default_rng(99)
n_genes_exch = 1000
y_exch = rng_exch.normal(size=n_genes_exch)
yc_exch = y_exch - y_exch.mean()
syy_exch = float(np.dot(yc_exch, yc_exch))

strata_signatures = [(5, 5), (20, 5), (10, 40), (50, 50), (3, 30)]
scores_by_stratum = {sig: [] for sig in strata_signatures}
for sig in strata_signatures:
    n_pos_e, n_neg_e = sig
    for _ in range(150):
        idx = rng_exch.choice(n_genes_exch, size=n_pos_e + n_neg_e, replace=False)
        pos_e, neg_e = idx[:n_pos_e].tolist(), idx[n_pos_e:].tolist()
        s, p, nu, b1 = tr._ulm_score(pos_e, neg_e, yc_exch, syy_exch, n_genes_exch)
        if s is not None:
            scores_by_stratum[sig].append(s)

pooled_sd = float(np.std(np.concatenate([np.array(v) for v in scores_by_stratum.values()])))
max_relative_gap = 0.0
for sig, vals in scores_by_stratum.items():
    stratum_sd = float(np.std(vals))
    gap = abs(stratum_sd - pooled_sd) / pooled_sd
    max_relative_gap = max(max_relative_gap, gap)
    print(f"    signature {sig}: n={len(vals)}  stratum_sd={stratum_sd:.3f}  pooled_sd={pooled_sd:.3f}")
check(max_relative_gap < 0.25,
      f"t-statistic spread is consistent (within 25%) across every (n_pos,n_neg) stratum "
      f"(max relative gap {max_relative_gap:.2%}) - the premise pooling all TFs' scores rests on")
check(abs(pooled_sd - 1.0) < 0.2, f"pooled sd is close to the theoretical t(n-2) sd of ~1.0 (got {pooled_sd:.3f})")


print("=== run_regulon_activity: calibration wiring is deterministic ===")
original_fetch = _patch_fetch(fake_regulons)
try:
    det_result_1 = tr.run_regulon_activity(resolved_up, ["TRRUST"], "9606", metric="log2fc", min_regulon_size=3, seed=7)
finally:
    tr.fetch_regulons = original_fetch
original_fetch = _patch_fetch(fake_regulons)
try:
    det_result_2 = tr.run_regulon_activity(resolved_up, ["TRRUST"], "9606", metric="log2fc", min_regulon_size=3, seed=7)
finally:
    tr.fetch_regulons = original_fetch
lam1 = det_result_1["tables"]["TRRUST"]["lambda_used"]
lam2 = det_result_2["tables"]["TRRUST"]["lambda_used"]
check(lam1 == lam2, f"identical run reproduces identical lambda_used (got {lam1} vs {lam2})")
p1_cal = det_result_1["tables"]["TRRUST"]["table"].iloc[0]["p-value (calibrated)"]
p2_cal = det_result_2["tables"]["TRRUST"]["table"].iloc[0]["p-value (calibrated)"]
check(p1_cal == p2_cal, f"identical run reproduces identical calibrated p-value (got {p1_cal} vs {p2_cal})")

print("=== run_regulon_activity: single-TF run correctly skips calibration (below min_tfs) ===")
skip_block = det_result_1["tables"]["TRRUST"]
check(skip_block["lambda_used"] == 1.0, "a 1-TF database run leaves lambda at 1.0 (below REGULON_MIN_TFS_FOR_CALIBRATION)")
check(skip_block["calibration_note"] is not None, "and explains why via calibration_note")


tr.REGULON_CACHE_DIR = _REAL_REGULON_CACHE_DIR
shutil.rmtree(_TEST_CACHE_DIR, ignore_errors=True)

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
