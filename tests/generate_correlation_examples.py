"""One-off generator: the three MATCHED-SAMPLE example tables for Multi-omics Correlation-based Integration.

    Data3/example_corr_transcriptomics.csv   genes (rows) x 12 samples (columns)        gene symbols
    Data3/example_corr_proteomics.csv        proteins (rows) x 12 samples (columns)      UniProt accessions
    Data3/example_corr_metabolomics.csv      12 samples (rows) x compounds + a Group     (a copy of the existing
                                             Data3/example_metabolomics_concentrations.csv - the 12 samples and
                                             their values are the anchor the other two are built around)

THESE ARE SYNTHETIC, on purpose, and the README says so: no public dataset with matched RNA, protein and metabolite
levels can be bundled, so the structure is PLANTED and known, which is what lets the tests assert what the method must
find rather than merely that it runs:

  - 15 TCA-cycle / glycolysis genes whose transcript and protein share a per-sample factor (a strong positive
    transcript-protein correlation, rho ~ 0.9);
  - 4 genes whose protein is the NEGATIVE of the transcript's factor (planted discordance, rho ~ -0.9);
  - CS, IDH3A and MDH2 are driven by the sample's citrate / isocitrate / malate level, so each is correlated with its own
    metabolite in all three omics (note: citrate, isocitrate and malate also differ between Control and Treated, so part of that
    correlation IS the group difference - exactly the caveat the README makes);
  - ~270 background genes whose RNA and protein are independent noise (no correlation, the null);
  - about 1% of the protein cells are missing.
The two gene tables use different layouts (samples in the columns; the metabolite file has samples in the rows with a
Group column) so the orientation detection is exercised by the example itself.

Seeded, so the files are reproducible. Needs the network once: the planted genes' UniProt accessions and the background
proteins' gene symbols come from the app's own resolvers (the accessions written here are exactly the ones the app
resolves back to the same symbols).
"""
import os
import shutil
import sys

import numpy as np
import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
import Genomics_module_retriever as gene_resolver  # noqa: E402
import Multiomics_module_retriever as mr  # noqa: E402

SEED = 20261007
N_BACKGROUND = 270
DATA = os.path.join(_ROOT, "Data3")


def main():
    rng = np.random.RandomState(SEED)
    src = os.path.join(DATA, "example_metabolomics_concentrations.csv")
    met = pd.read_csv(src)
    samples = list(met["Sample"])
    shutil.copyfile(src, os.path.join(DATA, "example_corr_metabolomics.csv"))
    n = len(samples)

    def z(col):
        v = met[col].to_numpy(dtype=float)
        return (v - v.mean()) / v.std(ddof=1)

    genes = list(pd.read_csv(os.path.join(DATA, "example_multiomics_genes.csv"))["Gene Symbol"])
    resolved = gene_resolver.resolve_gene_table(pd.DataFrame({"Gene Symbol": genes}), "gene_symbol", taxon_id="9606", verbose=False)
    acc = {r["Gene Symbol"]: r["UniProt ID"] for _, r in resolved.iterrows() if str(r["UniProt ID"]).strip()}
    planted = [g for g in genes if g in acc]
    concordant, discordant = planted[:15], planted[15:19]
    driven = {"CS": "Citrate", "IDH3A": "Isocitrate", "MDH2": "Malate"}

    rows_t, rows_p = {}, {}
    for g in concordant + discordant:
        f = z(driven[g]) if g in driven else rng.randn(n)
        sign = 1.0 if g in concordant else -1.0
        # noise is kept small on purpose: at 12 samples and thousands of pairs only a strong, clean relation survives the
        # multiple-testing correction (Spearman needs rho ~ 0.97 to reach p ~ 1e-7), and an example that shows nothing would
        # teach nothing
        t_noise, p_noise = (0.04, 0.08) if g in driven else (0.25, 0.35)
        rows_t[g] = 8.0 + 1.0 * (f + t_noise * rng.randn(n))
        rows_p[acc[g]] = 20.0 + 1.5 * (sign * f + p_noise * rng.randn(n))
    for g in planted[19:]:
        rows_t[g] = 8.0 + rng.randn(n)
        rows_p[acc[g]] = 20.0 + 1.5 * rng.randn(n)

    # background: proteins from the existing proteomics example, their gene symbols from the app's own resolver
    de = pd.read_csv(os.path.join(DATA, "example_proteomics_de_full.csv"))
    taken = set(acc.values())
    bg_acc = [a for a in de["Protein IDs"].astype(str) if a not in taken][: N_BACKGROUND * 2]
    bg = mr.resolve_matrix_features("proteomics", bg_acc, "9606")["gene_of"]
    kept, seen_sym = [], set(rows_t)
    for a in bg_acc:
        sym = bg.get(a)
        if sym and sym not in seen_sym:
            seen_sym.add(sym)
            kept.append((a, sym))
        if len(kept) == N_BACKGROUND:
            break
    for a, sym in kept:
        rows_t[sym] = 8.0 + rng.randn(n)
        rows_p[a] = 20.0 + 1.5 * rng.randn(n)

    tdf = pd.DataFrame(rows_t, index=samples).T.round(4)
    pdf = pd.DataFrame(rows_p, index=samples).T.round(4)
    miss = rng.choice(pdf.size, size=max(1, pdf.size // 100), replace=False)
    flat = pdf.to_numpy().copy()
    flat.ravel()[miss] = np.nan
    pdf = pd.DataFrame(flat, index=pdf.index, columns=pdf.columns)
    tdf.index.name, pdf.index.name = "Gene", "Protein IDs"
    tdf.to_csv(os.path.join(DATA, "example_corr_transcriptomics.csv"))
    pdf.to_csv(os.path.join(DATA, "example_corr_proteomics.csv"))
    print(f"transcriptomics: {tdf.shape[0]} genes x {tdf.shape[1]} samples; proteomics: {pdf.shape[0]} proteins x {pdf.shape[1]} samples "
          f"({int(pdf.isna().sum().sum())} missing); planted: {len(concordant)} concordant, {len(discordant)} discordant, "
          f"driven {sorted(driven)}; background {len(kept)}")


if __name__ == "__main__":
    main()
