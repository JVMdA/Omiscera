# -*- coding: utf-8 -*-
"""Evaluation Fix 3: production-scale FDR calibration + the permutation-count
sweep. A provenance script (prints tables, no PASS/FAIL bookkeeping - real
runtime, roughly an hour total), in the same spirit as probe_better_null_
estimator.py / calibration_regulon_activity.py: the measurements here justify
a design decision (whether Fix 1's skew-aware p-value is genuinely calibrated
at real scale, and what NETWORK_N_PERMUTATIONS_DEFAULT/CHOICES should be),
cited by name rather than re-run on every test pass.

PART A - TRUE NULL AT PRODUCTION SCALE. test_network_engine_core.py already
proves the FDR<0.05 call rate stays <=5% under a true null on 60-node random
graphs. Production is a real, structured ~16,201-protein STRING network with
heavy-tailed degree and genuine biological correlation between neighbours
(they share signal by construction) - that evidence does not transfer by
itself. This part repeats the identical experiment (seeds drawn uniformly at
random, no planted signal, so "observed" and "null" are from the same
distribution and every FDR<0.05 call is by definition a false positive) on
the REAL cached human graph, and reports the call rate under BOTH the
skew-aware (Fix 1) and the pre-fix normal-tail p-value - both derived from
the SAME permutation run's null moments/zscores, so no second (expensive)
permutation run is needed per repeat.

PART B - PERMUTATION-COUNT SWEEP. 200/500/1000/2000 permutations at a fixed
seed, on the bundled human DE example, reporting hits/recruited (both BH
families, see Fix 2)/measured null skewness at each - the evidence
NETWORK_N_PERMUTATIONS_DEFAULT/CHOICES (Transcriptomics_network.py) should be
set from, rather than the pre-evaluation guess of 500.
"""
import os
import sys
import time

_ROOT = r"c:\Users\cliente\OneDrive\Documentos\Doutorado_brasil\Resultados\app-metabolismo\Código\Omiscera_v06"
sys.path.insert(0, _ROOT)
import numpy as np
from scipy.special import ndtr

import Genomics_network as gnet
import Transcriptomics_module_retriever as retriever
import Transcriptomics_network as txnet

print("=" * 88)
print("PART A: true-null FDR calibration at PRODUCTION scale (real human STRING network)")
print("=" * 88)

graph = gnet.build_string_network("9606", confidence=700)
n_proteins = graph["n_proteins"]
print(f"real human network: {n_proteins} proteins\n")

N_REPEATS = 5
N_SEEDS = 550   # matches the real seed count this app typically sees on the bundled example
N_PERM = 500

skew_rates, normal_rates = [], []
for rep in range(N_REPEATS):
    rng = np.random.default_rng(2000 + rep)
    seed_indices = rng.choice(n_proteins, size=N_SEEDS, replace=False)
    seed_idx_score = {int(i): 1.0 for i in seed_indices}
    t0 = time.time()
    core = gnet._propagation_core(graph, seed_idx_score, n_permutations=N_PERM,
                                  seed_rng=int(rng.integers(1_000_000)))
    elapsed = time.time() - t0

    # Skew-aware (Fix 1, what is actually shipped).
    fdrs_skew = np.asarray(core["fdrs"])
    rate_skew = float((fdrs_skew < 0.05).mean())
    n_sig_skew = int((fdrs_skew < 0.05).sum())

    # Pre-fix normal-tail, recomputed from the SAME null moments/zscores -
    # not a second permutation run, just a different p-value formula applied
    # to the same observed data.
    pvals_normal = ndtr(-core["zscores"])
    fdrs_normal = np.asarray(gnet.genomics_ora.benjamini_hochberg(list(pvals_normal)))
    rate_normal = float((fdrs_normal < 0.05).mean())
    n_sig_normal = int((fdrs_normal < 0.05).sum())

    skew_rates.append(rate_skew)
    normal_rates.append(rate_normal)
    print(f"  rep {rep + 1}/{N_REPEATS} ({elapsed:.0f}s, n_gamma_tail={core['n_gamma_tail']}/{n_proteins}): "
          f"skew-aware {n_sig_skew}/{n_proteins} ({rate_skew:.4%}) significant  |  "
          f"pre-fix normal {n_sig_normal}/{n_proteins} ({rate_normal:.4%})")

print(f"\nSKEW-AWARE (Fix 1):  worst call rate {max(skew_rates):.4%} across {N_REPEATS} repeats - "
      f"{'CALIBRATED (<=5%)' if max(skew_rates) <= 0.05 else 'NOT CALIBRATED (>5%)'}")
print(f"PRE-FIX NORMAL:      worst call rate {max(normal_rates):.4%} across {N_REPEATS} repeats - "
      f"{'calibrated (<=5%)' if max(normal_rates) <= 0.05 else 'NOT CALIBRATED (>5%)'}")

print()
print("=" * 88)
print("PART B: permutation-count sweep (bundled human example, fixed seed=42)")
print("=" * 88)

de_df = retriever.load_de_table_from_bytes(
    open(os.path.join(_ROOT, "Data3", "example_transcriptomics_de_full.csv"), "rb").read())
resolved_df, _ = retriever.resolve_de_table(de_df, "9606", identifier_type="gene_symbol")

for n_perm in (200, 500, 1000, 2000):
    t0 = time.time()
    res = txnet.run_network_analysis(resolved_df, "9606", taxon_name="Homo sapiens",
                                     n_permutations=n_perm, seed=42, module_detection=False)
    elapsed = time.time() - t0
    prop = res["propagation"]
    skew = prop["null_diagnostics"]["null_skew_at_max"]
    n_gamma = prop["null_diagnostics"]["n_gamma_tail"]
    print(f"  n_permutations={n_perm:5d} ({elapsed:5.0f}s): hits={prop['n_hits']:4d}  "
          f"seed_hits={prop['n_hits_seed']:4d}  recruited={prop['n_recruited']:3d}  "
          f"recruited_all_family={prop['n_recruited_all_family']:3d}  "
          f"null_skew_at_max={skew:6.3f}  n_gamma_tail={n_gamma}/{prop['n_proteins']}")

print("\nDone.")
