"""LIVE check of KEGG in Proteomics ORA/Quantitative Enrichment - hits the
real KEGG REST API.

Not part of the offline suite (tests/test_proteomics_kegg.py covers the
logic with fixtures). This exists because the offline tests can only prove
the code does what the fixtures say; only a live run proves the accession ->
gene -> pathway join actually works against what KEGG serves, and that the
coverage numbers claimed in Proteomics_module.py's module docstring
(measured live before that code was written: 94.6% of a real 3,500-protein
human dataset maps to a unique KEGG gene) still hold.

Run manually:  python tests/live_proteomics_kegg.py
"""
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import pandas as pd

import Proteomics_module as pm

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


print("=== resolve_kegg_organism_code: reused unchanged, still resolves human/mouse ===")
# Full multi-organism coverage is already exercised by
# tests/live_transcriptomics_kegg.py against this SAME function (it is
# reused, not duplicated, for Proteomics) - this just confirms the import
# path Proteomics_module actually uses resolves correctly too.
for taxon_id, taxon_name, expected_code in [("9606", "Homo sapiens", "hsa"), ("10090", "Mus musculus", "mmu")]:
    res = pm.resolve_kegg_organism_code(taxon_id, taxon_name, force_refresh=True)
    check(res["organism_code"] == expected_code and res["verified"],
          f"{taxon_name} ({taxon_id}) -> {expected_code!r} via Proteomics_module's own import "
          f"(got {res['organism_code']!r}, verified={res['verified']})")

print("\n=== build_kegg_accession_sets: known human glycolysis/cell-cycle accessions ===")
GLYCOLYSIS_ACCESSIONS = {
    "P04406": "GAPDH", "P06733": "ENO1", "P14618": "PKM", "P04075": "ALDOA",
}
CELL_CYCLE_ACCESSIONS = {"P06493": "CDK1", "P14635": "CCNB1"}
measured = set(GLYCOLYSIS_ACCESSIONS) | set(CELL_CYCLE_ACCESSIONS) | {"NOTREAL00"}

result = pm.build_kegg_accession_sets("hsa", measured, force_refresh=True)
check(len(result["pathways"]) > 300, f"human has a full pathway set (got {len(result['pathways'])} pathways)")
check(result["fetch_failed"] is False, "clean fetch, no failure")
glyc = result["pathways"].get("hsa00010")
check(glyc is not None and glyc["name"] == "Glycolysis / Gluconeogenesis",
      f"hsa00010 keeps its real name (got {glyc['name'] if glyc else None})")
check(glyc is not None and set(GLYCOLYSIS_ACCESSIONS) <= glyc["genes"],
      f"all 4 real glycolysis accessions land in hsa00010's gene set "
      f"(got {glyc['genes'] & set(GLYCOLYSIS_ACCESSIONS) if glyc else None})")
cc = result["pathways"].get("hsa04110")
check(cc is not None and set(CELL_CYCLE_ACCESSIONS) <= cc["genes"],
      f"both real cell-cycle accessions land in hsa04110's gene set "
      f"(got {cc['genes'] & set(CELL_CYCLE_ACCESSIONS) if cc else None})")
check("NOTREAL00" not in (glyc["genes"] if glyc else set()) and result["n_mapped"] == 6,
      f"a fake accession maps to nothing and is excluded from n_mapped (got n_mapped={result['n_mapped']})")
check(result["n_measured"] == 7, f"n_measured counts the full measured set including the fake one (got {result['n_measured']})")

print("\n=== coverage on the app's own real example file (3,500 human proteins) ===")
example_path = os.path.join(_ROOT, "Data3", "example_proteomics_de_full.csv")
df = pd.read_csv(example_path)
id_col = df.columns[0]
accessions = {str(a).split(";")[0].strip() for a in df[id_col] if str(a).strip()}
check(len(accessions) > 3000, f"the example file has real scale (got {len(accessions)} accessions)")

cov = pm.build_kegg_accession_sets("hsa", accessions, force_refresh=True)
mapped_frac = cov["n_mapped"] / cov["n_measured"]
on_pathway_frac = cov["n_on_pathway"] / cov["n_measured"]
print(f"      n_measured={cov['n_measured']}, n_mapped={cov['n_mapped']} ({mapped_frac:.1%}), "
      f"n_on_pathway={cov['n_on_pathway']} ({on_pathway_frac:.1%})")
check(mapped_frac > 0.90,
      f"the documented ~94.6% accession->KEGG-gene coverage still holds on real data (got {mapped_frac:.1%})")
check(on_pathway_frac > 0.35,
      f"the documented ~44.1% pathway coverage still holds (got {on_pathway_frac:.1%})")
# The specific claim build_kegg_accession_sets' docstring makes about isoform
# accessions NOT double-counting a gene: every accession that lands in any
# pathway's set is still counted once per pathway (a Python set), and no
# accession here should silently vanish just for sharing a gene with another
# measured accession - spot check by re-deriving n_on_pathway from the
# pathway sets directly.
on_pathway_accs: set = set()
for info in cov["pathways"].values():
    on_pathway_accs |= info["genes"]
check(len(accessions & on_pathway_accs) == cov["n_on_pathway"],
      "n_on_pathway matches a direct recount from the returned pathway sets (no silent mismatch)")

print("\n=== end-to-end: run_ora + run_quantitative_enrichment against real human KEGG data ===")
rows = []
for acc in GLYCOLYSIS_ACCESSIONS:
    rows.append((acc, 3.0, 1e-8, 1e-7))
for acc in CELL_CYCLE_ACCESSIONS:
    rows.append((acc, -3.0, 1e-8, 1e-7))
# Padding so the universe isn't trivially small (also real human accessions,
# unrelated to either pathway, flat/non-significant).
OTHER = ["P60709", "P68371", "P62826", "P07437", "P11021", "P08238", "P30101",
         "P07900", "P0DMV8", "Q9Y230"]
for acc in OTHER:
    rows.append((acc, 0.05, 0.9, 0.95))

de = pd.DataFrame({
    "input_id": [r[0] for r in rows],
    "accession": [r[0] for r in rows],
    "gene_symbol": [GLYCOLYSIS_ACCESSIONS.get(r[0], CELL_CYCLE_ACCESSIONS.get(r[0], r[0])) for r in rows],
    "protein_name": [""] * len(rows),
    "log2FC": [r[1] for r in rows],
    "pvalue": [r[2] for r in rows],
    "padj": [r[3] for r in rows],
    "status": ["resolved"] * len(rows),
})

ora = pm.run_ora(
    de, databases=["KEGG"], taxon_id="9606", taxon_name="Homo sapiens",
    padj_cutoff=0.05, lfc_cutoff=1.0, min_term_size=2, max_term_size=500,
)
up = ora["tables"]["KEGG (up)"]["table"]
check(len(up) > 0, f"the up direction returned real KEGG terms (got {len(up)})")
check("hsa00010" in set(up["Term ID"].astype(str)),
      f"Glycolysis / Gluconeogenesis is enriched among the up-regulated proteins "
      f"(got {sorted(set(up['Term ID'].astype(str)))[:8]})")
down = ora["tables"]["KEGG (down)"]["table"]
check("hsa04110" in set(down["Term ID"].astype(str)),
      f"Cell cycle is enriched among the down-regulated proteins "
      f"(got {sorted(set(down['Term ID'].astype(str)))[:8]})")
overlap_col = up.loc[up["Term ID"] == "hsa00010", "Overlapping Proteins"]
check(len(overlap_col) and "P04406" in str(overlap_col.iloc[0]) and "GAPDH" in str(overlap_col.iloc[0]),
      f"reported members are the user's own accession+symbol, never a raw KEGG gene ID "
      f"(got {str(overlap_col.iloc[0])[:80] if len(overlap_col) else None!r})")

quant = pm.run_quantitative_enrichment(
    de, databases=["KEGG"], taxon_id="9606", taxon_name="Homo sapiens",
    min_term_size=2, max_term_size=500,
)
qtbl = quant["tables"]["KEGG"]["table"]
check(len(qtbl) > 0, f"Quantitative Enrichment tested real KEGG terms (got {len(qtbl)})")
effect_col = "Effect (mean difference)"
glyc_row = qtbl[qtbl["Term ID"] == "hsa00010"]
cc_row = qtbl[qtbl["Term ID"] == "hsa04110"]
check(len(glyc_row) > 0 and len(cc_row) > 0,
      f"both glycolysis and cell cycle were tested (got {len(glyc_row)}, {len(cc_row)} rows)")
if len(glyc_row) and len(cc_row):
    check(float(glyc_row[effect_col].iloc[0]) > float(cc_row[effect_col].iloc[0]),
          "glycolysis (up-regulated) scores higher than cell cycle (down-regulated) on the "
          f"Wilcoxon effect (got glycolysis={glyc_row[effect_col].iloc[0]}, "
          f"cell_cycle={cc_row[effect_col].iloc[0]})")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
