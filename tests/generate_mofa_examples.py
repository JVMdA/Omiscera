"""One-off generator: the three MATCHED-SAMPLE example tables for Multi-omics MOFA / Factor-based Integration.

    Data3/example_mofa_transcriptomics.csv   genes (rows) x 36 samples (columns)          gene symbols
    Data3/example_mofa_proteomics.csv        proteins (rows) x 36 samples (columns)        UniProt accessions (row i = the same gene as the RNA row i)
    Data3/example_mofa_metabolomics.csv      36 samples (rows) x 42 compounds + a Group    the compounds of example_metabolomics_concentrations.csv
    tests/fixtures/mofa_example_truth.csv    the TRUE factor values the data were built from (what the tests assert MOFA recovers)

THESE ARE SYNTHETIC, on purpose, and the README says so: no public dataset with matched RNA, protein and metabolite levels can be bundled, so the
structure is PLANTED and known:
  - FACTOR "treatment" (shared by all three omics): high in Treated, low in Control. Its genes are the members of Reactome's citric acid cycle
    (R-HSA-71403) and glycolysis (R-HSA-70171), in the transcripts AND the proteins; its metabolites are the TCA / glycolysis intermediates.
  - FACTOR "cholesterol" (transcripts and proteins ONLY, unrelated to the groups): the members of Reactome's cholesterol biosynthesis (R-HSA-191273).
  - FACTOR "lipids" (metabolites ONLY, unrelated to the groups): the lipid species of the panel (PC / PE / SM / TG / LPC / DG / ceramide, palmitic acid).
  - ~600 background genes that carry nothing but noise; about 4% of the protein cells are missing.
The pathways are real Reactome sets so the pathway enrichment of the weights has something true to find. Seeded, so the files are reproducible.
Needs the network once: the symbols of the chosen UniProt accessions come from the app's own resolver.
"""
import os
import sys

import numpy as np
import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
import Genomics_module as gm  # noqa: E402
import Multiomics_module_retriever as mr  # noqa: E402

SEED = 20261013
N_PER_GROUP = 18
N_BACKGROUND = 600
DATA = os.path.join(_ROOT, "Data3")
TCA_LIKE = ["Citrate", "Isocitrate", "2-Oxoglutarate", "Succinate", "Fumarate", "Malate", "Oxaloacetate", "cis-Aconitate", "Pyruvate", "Lactate",
            "Glucose 6-phosphate", "Fructose 6-phosphate", "Phosphoenolpyruvate", "3-Phosphoglycerate", "Dihydroxyacetone phosphate"]
LIPIDS = ["PC 34:1", "PC 36:2", "PC 32:0", "PE 34:1", "PE 36:2", "PE 38:4", "SM 34:1", "SM 36:1", "TG 52:2", "TG 54:3", "LPC 16:0", "LPC 18:0", "DG 36:2",
          "Cer d18:1/16:0", "Palmitic acid", "Carnitine"]


def main():
    rng = np.random.RandomState(SEED)
    uni = gm.build_reactome_universe("9606", taxon_name="Homo sapiens")["pathways"]
    f1_acc = sorted(set(uni["R-HSA-71403"]["genes"]) | set(uni["R-HSA-70171"]["genes"]))
    f2_acc = sorted(set(uni["R-HSA-191273"]["genes"]) - set(f1_acc))
    in_sets = set(f1_acc) | set(f2_acc)
    pool = sorted({a for info in uni.values() for a in info["genes"]} - in_sets)
    bg_acc = list(rng.choice(pool, size=N_BACKGROUND * 2, replace=False))
    accs = f1_acc + f2_acc + bg_acc
    print(f"treatment-factor genes {len(f1_acc)}, cholesterol-factor genes {len(f2_acc)}, candidate background {len(bg_acc)}")
    gene_of = mr.resolve_matrix_features("proteomics", accs, "9606")["gene_of"]
    rows, seen = [], set()
    for kind, lst in (("f1", f1_acc), ("f2", f2_acc), ("bg", bg_acc)):
        for a in lst:
            sym = gene_of.get(a)
            if sym and sym not in seen and (kind != "bg" or sum(1 for r in rows if r[0] == "bg") < N_BACKGROUND):
                seen.add(sym)
                rows.append((kind, a, sym))
    kinds = np.array([r[0] for r in rows])
    accs = [r[1] for r in rows]
    syms = [r[2] for r in rows]
    print(f"{len(rows)} genes: {(kinds == 'f1').sum()} treatment, {(kinds == 'f2').sum()} cholesterol, {(kinds == 'bg').sum()} background")

    samples = [f"Control_{i + 1}" for i in range(N_PER_GROUP)] + [f"Treated_{i + 1}" for i in range(N_PER_GROUP)]
    n = len(samples)
    group = np.array(["Control"] * N_PER_GROUP + ["Treated"] * N_PER_GROUP)
    z1 = (group == "Treated").astype(float) * 1.6 - 0.8 + 0.35 * rng.randn(n)
    z2 = rng.randn(n)
    z3 = rng.randn(n)
    z1, z2, z3 = ((z - z.mean()) / z.std() for z in (z1, z2, z3))

    def weights(mask, p, strength):
        w = np.where(mask, strength * (0.6 + 0.8 * rng.rand(p)), 0.0)
        return w * np.where(rng.rand(p) < 0.85, 1.0, -1.0)

    pg = len(rows)
    w1, w2 = weights(kinds == "f1", pg, 1.0), weights(kinds == "f2", pg, 1.0)
    rna = 8.0 + 1.0 * (np.outer(z1, w1) + np.outer(z2, w2)) + 0.8 * rng.randn(n, pg)
    w1p, w2p = w1 * (0.8 + 0.4 * rng.rand(pg)), w2 * (0.8 + 0.4 * rng.rand(pg))
    prot = 20.0 + 1.5 * (np.outer(z1, w1p) + np.outer(z2, w2p)) + 1.2 * rng.randn(n, pg)
    prot[rng.rand(n, pg) < 0.04] = np.nan

    comp = pd.read_csv(os.path.join(DATA, "example_metabolomics_concentrations.csv")).drop(columns=["Sample", "Group"]).columns.tolist()
    pm = len(comp)
    wm1 = np.array([(0.6 + 0.8 * rng.rand()) * (1 if rng.rand() < 0.85 else -1) if c in TCA_LIKE else 0.0 for c in comp])
    wm3 = np.array([(0.6 + 0.8 * rng.rand()) * (1 if rng.rand() < 0.85 else -1) if c in LIPIDS else 0.0 for c in comp])
    met = 30.0 + 5.0 * (np.outer(z1, wm1) + np.outer(z3, wm3)) + 4.0 * rng.randn(n, pm)
    met = np.round(np.clip(met, 0.5, None), 3)

    os.makedirs(DATA, exist_ok=True)
    tdf = pd.DataFrame(np.round(rna, 4).T, index=syms, columns=samples)
    tdf.index.name = "Gene"
    pdf = pd.DataFrame(np.round(prot, 4).T, index=accs, columns=samples)
    pdf.index.name = "Protein IDs"
    tdf.to_csv(os.path.join(DATA, "example_mofa_transcriptomics.csv"))
    pdf.to_csv(os.path.join(DATA, "example_mofa_proteomics.csv"))
    mdf = pd.DataFrame(met, columns=comp)
    mdf.insert(0, "Group", group)
    mdf.insert(0, "Sample", samples)
    mdf.to_csv(os.path.join(DATA, "example_mofa_metabolomics.csv"), index=False)
    os.makedirs(os.path.join(_ROOT, "tests", "fixtures"), exist_ok=True)
    pd.DataFrame({"Sample": samples, "Group": group, "treatment": z1, "cholesterol": z2, "lipids": z3}).to_csv(
        os.path.join(_ROOT, "tests", "fixtures", "mofa_example_truth.csv"), index=False)
    print(f"transcriptomics {tdf.shape}, proteomics {pdf.shape} ({int(np.isnan(prot).sum())} missing), metabolomics {mdf.shape}")


if __name__ == "__main__":
    main()
