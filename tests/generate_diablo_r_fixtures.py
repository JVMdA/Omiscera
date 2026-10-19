"""Regenerates tests/fixtures/diablo_reference/{block_*.csv, group.csv, train/test_samples.csv}: the SMALL, matched
sample x feature tables (36 samples, top 60 most-variable transcripts / proteins, the 42 metabolites, prepared
EXACTLY as Multiomics_diablo.prepare_blocks would - dropped constant / too-missing features, kept top-variance) that
tests/r_reference/diablo_reference.R feeds to a real mixOmics block.splsda, together with a fixed, stratified
70/30 train/test split. Reuses the existing Data3/example_mofa_*.csv tables (36 samples, Control/Treated,
missing protein values) - no new example data, so Multi-omics MOFA and DIABLO look at the same uploads.

After running this (offline, no network needed - the tables are already on disk), re-run
tests/r_reference/diablo_reference.R (needs a local R + mixOmics - see tests/r_reference/README.txt) to refresh
tests/fixtures/diablo_reference/out/*.csv, the fixed mixOmics output tests/test_multiomics_diablo.py checks the
engine against. Not itself a test.
"""
import os
import sys

import numpy as np
import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
import Multiomics_correlation as mc  # noqa: E402

OUT = os.path.join(_ROOT, "tests", "fixtures", "diablo_reference")
SEED = 777
TOP_VAR = 60


def main():
    os.makedirs(OUT, exist_ok=True)
    mats = {}
    for o, f in (("transcriptomics", "example_mofa_transcriptomics.csv"), ("proteomics", "example_mofa_proteomics.csv"),
                 ("metabolomics", "example_mofa_metabolomics.csv")):
        mats[o] = mc.parse_matrix(open(os.path.join(_ROOT, "Data3", f), "rb").read())["df"]
    al = mc.align_samples(mats)
    samples = al["samples"]
    group = pd.read_csv(os.path.join(_ROOT, "Data3", "example_mofa_metabolomics.csv")).set_index("Sample").loc[samples, "Group"]

    rng = np.random.RandomState(SEED)
    for o in mats:
        out, _ = mc.prepare_matrix(al["mats"][o], log2=False, max_missing=0.5, top_variance=TOP_VAR)
        out.to_csv(os.path.join(OUT, f"block_{o}.csv"))
        print(o, out.shape, "missing cells:", int(out.isna().sum().sum()))
    group.to_csv(os.path.join(OUT, "group.csv"))

    idx = np.arange(len(samples))
    train_idx, test_idx = [], []
    for g in sorted(group.unique()):
        gi = idx[group.to_numpy() == g]
        rng.shuffle(gi)
        k = max(2, round(len(gi) * 0.7))
        train_idx += list(gi[:k])
        test_idx += list(gi[k:])
    train_idx, test_idx = sorted(train_idx), sorted(test_idx)
    pd.DataFrame({"sample": [samples[i] for i in train_idx]}).to_csv(os.path.join(OUT, "train_samples.csv"), index=False)
    pd.DataFrame({"sample": [samples[i] for i in test_idx]}).to_csv(os.path.join(OUT, "test_samples.csv"), index=False)
    print("train", len(train_idx), "test", len(test_idx))


if __name__ == "__main__":
    main()
