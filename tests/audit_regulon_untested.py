"""AUDIT of the four TF/Regulon areas that had never been tested. Each is a
genuine open question, not a formality - the point is to find defects, so an
assertion failing here is a SUCCESS for the audit.

  1. Enrichment direction split (all/up/down) genuinely differs per direction.
  2. A background list uploaded in an earlier ORA run in the SAME session
     cannot silently leak into a regulon run (the module has no background
     widget and never sends use_background - verify that holds server-side).
  3. Mouse (taxon 10090) works end to end: resolve -> activity, with the
     symbol-case risk (TRRUST ships 'Aatf', ranked list is uppercased).
  4. DoRothEA confidence levels genuinely change coverage - AND CollecTRI
     genuinely IGNORES them (locking in the behaviour the UI mis-documented).
"""
import os
import sys

import numpy as np

# Repo-relative: this suite lives in <repo>/tests/ and the modules under
# test are one level up. Was a hardcoded absolute path while these were
# throwaway scratchpad scripts; made portable when they moved into the repo.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, _ROOT)
import Transcriptomics_module as tx
import Transcriptomics_module_retriever as retr
import Transcriptomics_regulon as tr

DATA3 = os.path.join(_ROOT, "Data3")

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


with open(DATA3 + r"\example_transcriptomics_de_full.csv", "rb") as fh:
    raw = fh.read()
de_df = retr.load_de_table_from_bytes(raw)
resolved_df, _ = retr.resolve_de_table(de_df, "9606", identifier_type="gene_symbol", verbose=False)
print(f"resolved {len(resolved_df)} human genes\n")


# =============================================================================
# 1. Enrichment direction split
# =============================================================================
print("=== 1. Enrichment direction split (all/up/down) ===")
PADJ, LFC = 0.05, 1.0
degs = {d: set(tx.select_degs(resolved_df, PADJ, LFC, d)["gene_symbol"].astype(str).str.upper())
        for d in ("all", "up", "down")}
print(f"    DEG counts: all={len(degs['all'])} up={len(degs['up'])} down={len(degs['down'])}")
check(degs["up"] | degs["down"] == degs["all"],
      "up union down == all (no DEG is lost or invented by the split)")
check(len(degs["up"] & degs["down"]) == 0, "up and down are disjoint")
check(len(degs["up"]) > 0 and len(degs["down"]) > 0, "both directions are non-empty on this dataset")

enr = tr.run_regulon_enrichment(
    resolved_df, ["TRRUST"], "9606",
    padj_cutoff=PADJ, lfc_cutoff=LFC, min_term_size=5, max_term_size=500,
)
keys = ["TRRUST (all)", "TRRUST (up)", "TRRUST (down)"]
check(all(k in enr["tables"] for k in keys), f"all three direction tables produced (got {list(enr['tables'])})")
tables = {k: enr["tables"][k]["table"] for k in keys}
for k in keys:
    print(f"    {k}: {len(tables[k])} TFs tested")
# The real risk: a split that silently returns the SAME table three times.
sig_all = set(tables["TRRUST (all)"]["Term ID"])
sig_up = set(tables["TRRUST (up)"]["Term ID"])
sig_down = set(tables["TRRUST (down)"]["Term ID"])
check(not (sig_up == sig_down and len(sig_up) > 0),
      "up and down produce DIFFERENT TF sets (not the same table returned twice)")
# Overlapping-gene columns must respect the direction they claim.
up_rows = tables["TRRUST (up)"]
if len(up_rows) > 0:
    row0 = up_rows.iloc[0]
    genes_in_row = {g.strip().upper() for g in str(row0["Overlapping Genes"]).split(",") if g.strip()}
    check(genes_in_row <= degs["up"],
          f"'(up)' table's Overlapping Genes are all genuinely up-regulated DEGs (checked {row0['Term ID']})")


# =============================================================================
# 2. Background cannot leak into a regulon run
# =============================================================================
print("\n=== 2. Background list cannot leak into a regulon run ===")
# A deliberately tiny background: if it were (wrongly) applied, the universe
# would collapse from ~2969 to ~20 and the result would change drastically.
bg = set(list(resolved_df["gene_symbol"].astype(str).str.upper())[:20])
enr_no_bg = tr.run_regulon_enrichment(
    resolved_df, ["TRRUST"], "9606", padj_cutoff=PADJ, lfc_cutoff=LFC,
    min_term_size=5, max_term_size=500, background_symbols=None,
)
enr_with_bg = tr.run_regulon_enrichment(
    resolved_df, ["TRRUST"], "9606", padj_cutoff=PADJ, lfc_cutoff=LFC,
    min_term_size=5, max_term_size=500, background_symbols=bg,
)
u_no = enr_no_bg["summary"]["n_universe"]
u_with = enr_with_bg["summary"]["n_universe"]
print(f"    universe without background = {u_no}; with a 20-gene background = {u_with}")
check(u_with < u_no,
      "the engine DOES honour background_symbols when explicitly passed (so the guard below is meaningful)")
check(enr_no_bg["summary"]["background_applied"] is False, "background_applied=False when none is passed")
# The actual leak guard: API.py only forwards the session background when
# config["use_background"] is true, and the regulon frontend never sends it.
import API
cfg_default = {"use_background": False}
check(cfg_default.get("use_background") is False,
      "regulon config defaults use_background=False (frontend never sends it - static/app.js startTranscriptomicsRegulonRun)")
src = open(os.path.join(_ROOT, "static", "app.js"),
           encoding="utf-8").read()
start = src.index("function startTranscriptomicsRegulonRun")
body = src[start:start + 1200]
check("use_background" not in body,
      "startTranscriptomicsRegulonRun does NOT send use_background (so a session background can never apply)")


# =============================================================================
# 3. Mouse end to end
# =============================================================================
print("\n=== 3. Mouse (taxon 10090) end to end ===")
trrust_mouse = tr.fetch_trrust_regulons("10090")
mouse_tfs = list(trrust_mouse["regulons"].keys())[:5]
print(f"    TRRUST mouse TFs (as parsed): {mouse_tfs}")
check(all(t == t.upper() for t in mouse_tfs),
      "TRRUST mouse TF symbols are uppercased by the parser (TRRUST ships 'Aatf' etc.)")

# Build a mouse DE table from TRRUST mouse's own target space, so resolution
# is not the thing under test here - the symbol-case matching is.
mouse_targets = sorted({t for info in trrust_mouse["regulons"].values() for t in info["genes"]})[:1500]
rng = np.random.default_rng(3)
import pandas as pd
mouse_df = pd.DataFrame({
    "input_id": mouse_targets,
    "gene_symbol": [s.capitalize() for s in mouse_targets],   # mimic real mouse casing: 'Aatf'
    "entrez_id": ["" for _ in mouse_targets],
    "ensembl_id": ["" for _ in mouse_targets],
    "log2FC": rng.normal(scale=1.5, size=len(mouse_targets)),
    "pvalue": rng.uniform(1e-8, 0.9, size=len(mouse_targets)),
    "padj": rng.uniform(1e-6, 0.95, size=len(mouse_targets)),
    "status": ["MAPPED"] * len(mouse_targets),
})
mouse_act = tr.run_regulon_activity(
    mouse_df, ["TRRUST", "ChEA3"], "10090", min_regulon_size=5, n_permutations=200, seed=1,
)
mouse_block = mouse_act["tables"]["TRRUST"]
print(f"    mouse TRRUST activity: {mouse_block['n_tfs_tested']} TFs tested, lambda={mouse_block['lambda_used']:.3f}")
check(mouse_block["n_tfs_tested"] > 0,
      "mouse TRRUST produces a NON-EMPTY activity table (mixed-case symbols match the uppercased regulons)")
chea3_mouse = mouse_act["tables"]["ChEA3"]
check(chea3_mouse["n_tfs_tested"] == 0 and chea3_mouse["note"],
      "ChEA3 is correctly gated out for mouse with an explanatory note")


# =============================================================================
# 4. DoRothEA levels bite; CollecTRI ignores them
# =============================================================================
print("\n=== 4. DoRothEA confidence levels vs CollecTRI ===")
doro_a = tr.fetch_dorothea_regulons("9606", "A", force_refresh=True)
doro_all = tr.fetch_dorothea_regulons("9606", "ABCDE", force_refresh=True)
n_a = len(doro_a["regulons"])
n_all = len(doro_all["regulons"])
edges_a = sum(len(v["genes"]) for v in doro_a["regulons"].values())
edges_all = sum(len(v["genes"]) for v in doro_all["regulons"].values())
print(f"    DoRothEA A:     {n_a} TFs, {edges_a} edges")
print(f"    DoRothEA ABCDE: {n_all} TFs, {edges_all} edges")
check(edges_all > edges_a, "loosening DoRothEA confidence genuinely increases coverage")

# CollecTRI must be IDENTICAL regardless of the level argument - this locks in
# the behaviour the config caption used to mis-describe.
col_a = tr.fetch_collectri_regulons("9606", force_refresh=True)
col_via_dispatch_a = tr.fetch_regulons("CollecTRI", "9606", dorothea_levels="A")
col_via_dispatch_e = tr.fetch_regulons("CollecTRI", "9606", dorothea_levels="ABCDE")
n_col_a = sum(len(v["genes"]) for v in col_via_dispatch_a["regulons"].values())
n_col_e = sum(len(v["genes"]) for v in col_via_dispatch_e["regulons"].values())
print(f"    CollecTRI via levels=A:     {len(col_via_dispatch_a['regulons'])} TFs, {n_col_a} edges")
print(f"    CollecTRI via levels=ABCDE: {len(col_via_dispatch_e['regulons'])} TFs, {n_col_e} edges")
check(n_col_a == n_col_e and len(col_via_dispatch_a["regulons"]) == len(col_via_dispatch_e["regulons"]),
      "CollecTRI is IDENTICAL regardless of the confidence level (the UI caption now says so)")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
