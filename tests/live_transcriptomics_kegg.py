"""LIVE check of KEGG in Transcriptomics ORA/GSEA — hits the real KEGG REST API.

Not part of the offline suite (tests/test_transcriptomics_kegg.py covers the
logic with fixtures). This exists because the offline tests can only prove the
code does what the fixtures say; only a live run proves the fixtures match
what KEGG actually serves, and that the two bridges KEGG needs here —
organism (taxon -> KEGG organism code) and identifier (KEGG gene -> Entrez ->
this run's gene symbols) — line up against real data for a model organism AND
for a non-model one, which is the whole reason KEGG is worth offering.

Run manually:  python tests/live_transcriptomics_kegg.py
"""
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import pandas as pd

import Transcriptomics_module as tm

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


print("=== resolve_kegg_organism_code against the real KEGG genome list ===")
EXPECTED = [
    ("9606", "Homo sapiens", "hsa"),
    ("10090", "Mus musculus", "mmu"),
    ("10116", "Rattus norvegicus", "rno"),
    ("7955", "Danio rerio", "dre"),
    ("7227", "Drosophila melanogaster", "dme"),
    ("6239", "Caenorhabditis elegans", "cel"),
    ("3702", "Arabidopsis thaliana", "ath"),
    # Non-model: the reason KEGG is offered here at all.
    ("3827", "Cicer arietinum", "cam"),
    # The app's organism name carries a strain suffix KEGG's genome name does
    # not — resolution must fall back to genus+species and still verify.
    ("559292", "Saccharomyces cerevisiae S288C", "sce"),
]
for taxon_id, taxon_name, expected_code in EXPECTED:
    res = tm.resolve_kegg_organism_code(taxon_id, taxon_name, force_refresh=True)
    check(res["organism_code"] == expected_code and res["verified"],
          f"{taxon_name} ({taxon_id}) -> {expected_code!r}, taxon-verified "
          f"(got {res['organism_code']!r}, verified={res['verified']})")

print("\n=== an organism KEGG cannot resolve unambiguously ===")
res = tm.resolve_kegg_organism_code("1423", "Bacillus subtilis", force_refresh=True)
check(res["organism_code"] == "" and res["note"],
      "SENSITIVITY: species-level taxon 1423 matches many strain-level KEGG genomes and "
      f"none carries taxon 1423 -> refuses rather than guessing (got {res['organism_code']!r})")
print("      note: " + (res["note"] or "")[:160])
check(len(res["candidates"]) >= 10,
      f"the real KEGG genome list carries this many B. subtilis strains as structured "
      f"candidates too, not just in the note text (got {len(res['candidates'])})")
codes = {c["organism_code"] for c in res["candidates"]}
check("bsu" in codes,
      f"the canonical strain 'bsu' (B. subtilis 168) is among the real candidates (got "
      f"{sorted(codes)[:8]})")

print("\n=== confirming a real B. subtilis candidate against live KEGG ===")
confirmed = tm.resolve_kegg_organism_code(
    "1423", "Bacillus subtilis", force_refresh=True, confirmed_organism_code="bsu",
)
check(confirmed["organism_code"] == "bsu" and confirmed["verified"] is False,
      f"a real candidate is accepted as a manual confirmation, still marked unverified "
      f"(got {confirmed['organism_code']!r}, verified={confirmed.get('verified')})")
bsu_sets = tm.fetch_kegg_sets("bsu", force_refresh=True)
check(len(bsu_sets["pathways"]) > 50,
      f"the confirmed organism code fetches a real KEGG pathway set, exactly like any "
      f"automatically-resolved organism would (got {len(bsu_sets['pathways'])} pathways)")

print("\n=== confirming a code KEGG never actually matched by name is rejected ===")
rejected = tm.resolve_kegg_organism_code(
    "1423", "Bacillus subtilis", force_refresh=True, confirmed_organism_code="hsa",
)
check(rejected["organism_code"] == "",
      f"'hsa' is not one of B. subtilis's real name-matched candidates, so it is ignored "
      f"against the live candidate list too, not just the offline fixture "
      f"(got {rejected['organism_code']!r})")

print("\n=== fetch_kegg_sets: human ===")
hsa = tm.fetch_kegg_sets("hsa", force_refresh=True)
check(len(hsa["pathways"]) > 300,
      f"human has a full pathway set (got {len(hsa['pathways'])} pathways)")
check(hsa["caveat"] is None and not hsa["fetch_failed"], "clean fetch, no caveat")
bridged_frac = len(hsa["pathways"]) / max(hsa["n_pathways_before_bridge"], 1)
check(bridged_frac > 0.95,
      f"nearly every KEGG pathway survives the KEGG-gene -> Entrez bridge "
      f"({len(hsa['pathways'])}/{hsa['n_pathways_before_bridge']} = {bridged_frac:.1%})")
glyc = hsa["pathways"].get("hsa00010")
check(glyc is not None and glyc["name"] == "Glycolysis / Gluconeogenesis",
      f"hsa00010 keeps its real name, organism suffix stripped (got "
      f"{glyc['name'] if glyc else None})")
check(glyc is not None and "2597" in glyc["genes"],
      "and contains GAPDH's Entrez ID (2597) — the KEGG gene -> Entrez conversion is real, "
      f"not an assumed identity (glycolysis has {len(glyc['genes']) if glyc else 0} genes)")
check(all(g.isdigit() for g in glyc["genes"]),
      "every member is a bare numeric Entrez ID, with no 'hsa:' prefix left over")

print("\n=== fetch_kegg_sets: chickpea (non-model) ===")
cam = tm.fetch_kegg_sets("cam", force_refresh=True)
check(len(cam["pathways"]) > 100,
      f"chickpea has real KEGG pathway coverage (got {len(cam['pathways'])} pathways) — the "
      "case GO/Reactome/MSigDB serve poorly and KEGG is here for")
cam_genes = {g for info in cam["pathways"].values() for g in info["genes"]}
check(all(g.isdigit() for g in cam_genes),
      f"and its {len(cam_genes)} pathway genes are all bare Entrez IDs too")

print("\n=== end-to-end: run_ora + run_gsea against real human KEGG data ===")
# A small, deliberately glycolysis-heavy DE table: the point is not the
# biology but that real KEGG Entrez IDs meet real gene symbols through the
# resolved table's own bridge.
GLYCOLYSIS = [("GAPDH", "2597"), ("ENO1", "2023"), ("PKM", "5315"), ("ALDOA", "226"),
              ("PGK1", "5230"), ("TPI1", "7167"), ("HK1", "3098"), ("PFKL", "5211")]
CELL_CYCLE = [("CDK1", "983"), ("CCNB1", "891"), ("CDC20", "991"), ("BUB1", "699"),
              ("PLK1", "5347"), ("MAD2L1", "4085")]
OTHER = [("ACTB", "60"), ("TUBB", "203068"), ("RPL13A", "23521"), ("B2M", "567"),
         ("HPRT1", "3251"), ("SDHA", "6389"), ("YWHAZ", "7534"), ("PPIA", "5478")]

rows = []
for sym, entrez in GLYCOLYSIS:
    rows.append((sym, entrez, 3.0, 1e-8, 1e-7))
for sym, entrez in CELL_CYCLE:
    rows.append((sym, entrez, -3.0, 1e-8, 1e-7))
for sym, entrez in OTHER:
    rows.append((sym, entrez, 0.05, 0.9, 0.95))

de = pd.DataFrame({
    "input_id": [r[0] for r in rows],
    "gene_symbol": [r[0] for r in rows],
    "entrez_id": [r[1] for r in rows],
    "ensembl_id": [""] * len(rows),
    "log2FC": [r[2] for r in rows],
    "pvalue": [r[3] for r in rows],
    "padj": [r[4] for r in rows],
    "status": ["resolved"] * len(rows),
})

ora = tm.run_ora(
    de, databases=["KEGG"], taxon_id="9606", taxon_name="Homo sapiens",
    padj_cutoff=0.05, lfc_cutoff=1.0, min_term_size=3, max_term_size=500,
    use_length_correction=False,
)
up = ora["tables"]["KEGG (up)"]["table"]
check(len(up) > 0, f"the up direction returned KEGG terms (got {len(up)})")
top_up = set(up.head(5)["Term ID"].astype(str))
check("hsa00010" in top_up,
      f"Glycolysis / Gluconeogenesis (hsa00010) is among the top up-regulated terms "
      f"(top 5: {sorted(top_up)})")
down = ora["tables"]["KEGG (down)"]["table"]
check("hsa04110" in set(down.head(5)["Term ID"].astype(str)),
      f"Cell cycle (hsa04110) is among the top down-regulated terms "
      f"(top 5: {sorted(set(down.head(5)['Term ID'].astype(str)))})")
name_col = up.loc[up["Term ID"] == "hsa00010", "Term Name"]
check(len(name_col) and " - Homo sapiens" not in str(name_col.iloc[0]),
      f"the displayed term name carries no leftover organism suffix "
      f"(got {name_col.iloc[0] if len(name_col) else None!r})")
genes_col = up.loc[up["Term ID"] == "hsa00010", "Overlapping Genes"]
check(len(genes_col) and "GAPDH" in str(genes_col.iloc[0]),
      "and the reported genes are the user's own symbols, never KEGG gene IDs or Entrez IDs "
      f"(got {str(genes_col.iloc[0])[:70] if len(genes_col) else None})")

gsea = tm.run_gsea(
    de, databases=["KEGG"], taxon_id="9606", taxon_name="Homo sapiens",
    min_term_size=3, max_term_size=500, n_permutations=1000, seed=42,
)
gtable = gsea["tables"]["KEGG"]["table"]
check(len(gtable) > 0, f"GSEA tested real KEGG terms (got {len(gtable)})")
nes = dict(zip(gtable["Term ID"].astype(str), gtable["NES"]))
check(nes.get("hsa00010", 0) > 0,
      f"glycolysis gets a positive NES against real KEGG membership (got {nes.get('hsa00010')})")
check(nes.get("hsa04110", 0) < 0,
      f"cell cycle a negative one (got {nes.get('hsa04110')})")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
