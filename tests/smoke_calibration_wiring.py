"""Smoke-check the calibration wiring end-to-end through run_regulon_activity
on real data: columns present, FDR sorted by calibrated, lambda/note/center
surfaced on the block, and calibrated FDR is >= parametric FDR in aggregate
(calibration should only ever be as strict or stricter)."""
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

with open(os.path.join(_ROOT, "Data3", "example_transcriptomics_de_full.csv"), "rb") as fh:
    raw = fh.read()
de_df = retr.load_de_table_from_bytes(raw)
resolved_df, _ = retr.resolve_de_table(de_df, "9606", identifier_type="gene_symbol", verbose=False)

result = tr.run_regulon_activity(
    resolved_df, ["TRRUST", "DoRothEA", "CollecTRI", "ChEA3"], "9606",
    metric=tx.GSEA_RANKING_METRIC_DEFAULT, min_regulon_size=5, n_permutations=500, seed=42,
)

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


for db, block in result["tables"].items():
    print(f"\n--- {db} ---")
    print(f"  lambda_used={block['lambda_used']:.4f}  null_center={block['null_center']:.4f}  note={block['calibration_note']}")
    df = block["table"]
    if df.empty:
        check(block["calibration_note"] is not None or db == "ChEA3",
              f"{db}: empty table has an explanatory note or is the expected unsigned skip")
        continue
    check(list(df.columns) == tr.ACTIVITY_RESULT_COLUMNS, f"{db}: columns match ACTIVITY_RESULT_COLUMNS")
    # Sorted by calibrated FDR ascending.
    fdrs = df["FDR (calibrated)"].to_numpy(dtype=float)
    check(np.all(np.diff(fdrs) >= -1e-12), f"{db}: table sorted ascending by FDR (calibrated)")
    # z = (score - center) / lambda, recomputed independently.
    recomputed_z = (df["Activity Score"].to_numpy(dtype=float) - block["null_center"]) / block["lambda_used"]
    check(np.allclose(recomputed_z, df["Calibrated Score (z)"].to_numpy(dtype=float), atol=1e-9),
          f"{db}: Calibrated Score (z) matches (score-center)/lambda recomputed independently")
    # Calibration only ever widens or matches the parametric p-value's rank
    # order of confidence in aggregate - not a per-row guarantee (z-test vs
    # t-test differ), but median calibrated p should be >= median parametric
    # p whenever lambda > 1.
    if block["lambda_used"] > 1.0 + 1e-9:
        med_param = float(np.median(df["p-value (parametric)"]))
        med_cal = float(np.median(df["p-value (calibrated)"]))
        check(med_cal >= med_param - 1e-9,
              f"{db}: lambda={block['lambda_used']:.3f}>1, median calibrated p ({med_cal:.4g}) >= median parametric p ({med_param:.4g})")
    n_sig_param = int((df["FDR (parametric)"] < 0.05).sum())
    n_sig_cal = int((df["FDR (calibrated)"] < 0.05).sum())
    print(f"  n_tfs={len(df)}  sig@FDR<0.05 parametric={n_sig_param}  calibrated={n_sig_cal}")
    if block["lambda_used"] > 1.0 + 1e-9:
        check(n_sig_cal <= n_sig_param,
              f"{db}: calibration (lambda>1) does not INCREASE the count of significant TFs ({n_sig_cal} <= {n_sig_param})")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
