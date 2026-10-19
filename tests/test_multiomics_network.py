# -*- coding: utf-8 -*-
"""Offline tests of Multiomics_network.py (Multi-omics Network-based Integration).

What is checked, each against an INDEPENDENT reference or a hand computation:
  - the equal-weight-per-omics rule (hand-computed weights; ONE omics gives exactly 1.0);
  - metabolite -> enzyme: a brute-force join of two small files, the skeleton bridge, the hub guard at its boundary;
  - the SINGLE MOST IMPORTANT CHECK: with the same seeds the combined network equals the existing Proteomics network tile
    exactly (same propagation, same null, same FDR) - so nothing here can drift from the numerics the app already trusts;
  - several omics end to end on a 3-community STRING-like graph: who seeds each protein, per-omics module columns;
  - the pathway-overlap network against a brute-force Jaccard (threshold boundary, compounds kept apart from genes,
    components, top-N), and a layout in which no two nodes overlap.
"""
import gzip
import math
import os
import sys
import tempfile

import numpy as np
import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import Genomics_module as gm  # noqa: E402
import Genomics_network as gnet  # noqa: E402
import Multiomics_network as mn  # noqa: E402
import Proteomics_network as pnet  # noqa: E402

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def raises(fn, *a, **k):
    try:
        fn(*a, **k)
    except Exception as exc:
        return type(exc).__name__, str(exc)
    return None, ""


# =============================================================================
print("=== equal weight per omics ===")
# =============================================================================
one = mn.combine_seed_weights({"proteomics": {"a", "b", "c"}})
check(all(abs(w - 1.0) < 1e-12 for w in one["weights"].values()) and one["omics"] == ["proteomics"],
      "ONE omics: every seed weighs exactly 1.0 (the single-omics tile)")
two = mn.combine_seed_weights({"transcriptomics": {"a", "b"}, "proteomics": {"b", "c", "d", "e"}})
# n = 5, K = 2: w(p) = 5 * 1/2 * sum_o [p in S_o]/|S_o|
exp = {"a": 5 * 0.5 * (1 / 2), "b": 5 * 0.5 * (1 / 2 + 1 / 4), "c": 5 * 0.5 * (1 / 4), "d": 5 * 0.5 * (1 / 4), "e": 5 * 0.5 * (1 / 4)}
check(all(abs(two["weights"][p] - exp[p]) < 1e-12 for p in exp) and two["n_seeds"] == 5,
      f"hand-computed weights for two overlapping seed sets ({ {p: round(v, 4) for p, v in two['weights'].items()} })")
check(abs(sum(two["weights"].values()) - 5) < 1e-12, "the weights sum to the number of distinct seeds (mean weight 1.0)")
check(two["seed_omics"]["b"] == ["transcriptomics", "proteomics"] and two["seed_omics"]["a"] == ["transcriptomics"],
      "each protein records which omics seeded it, in the fixed omics order")
check(abs(two["weights"]["b"] - two["weights"]["a"] - 5 * 0.5 * (1 / 4)) < 1e-12 and two["weights"]["b"] > two["weights"]["a"],
      "SENSITIVITY: a protein hit by BOTH omics outweighs one hit by one")
big_small = mn.combine_seed_weights({"transcriptomics": {f"t{i}" for i in range(3000)}, "metabolomics": {"m1", "m2", "m3", "m4"}})
share_t = sum(w for p, w in big_small["weights"].items() if p.startswith("t"))
share_m = sum(w for p, w in big_small["weights"].items() if p.startswith("m"))
check(abs(share_t - share_m) < 1e-9, f"SENSITIVITY: 3,000 transcript seeds and 4 metabolite seeds carry EQUAL total weight ({share_t:.3f} vs {share_m:.3f})")
three = mn.combine_seed_weights({"genomics": {"x"}, "transcriptomics": set(), "proteomics": {"x", "y"}})
check(three["omics"] == ["genomics", "proteomics"] and abs(three["weights"]["x"] - 2 * 0.5 * (1 + 0.5)) < 1e-12,
      "an omics with no seeds is ignored, not counted in K")
check(mn.combine_seed_weights({})["n_seeds"] == 0, "no seeds at all gives an empty result, not an error")

# =============================================================================
print("=== metabolite -> enzyme ===")
# =============================================================================
CHEBI_RX = {"1": ["R1", "R2"], "2": ["R2", "R3", "R4"], "3": ["R5"], "4": [f"H{i}" for i in range(61)], "5": ["R9"], "6": [f"G{i}" for i in range(60)]}
RX_PROT = {"R1": ["P1", "P2"], "R2": ["P2"], "R3": ["P3"], "R4": ["P1"], "R5": ["P5"], "R9": []}
RX_PROT.update({f"H{i}": [f"HP{i}"] for i in range(61)})
RX_PROT.update({f"G{i}": [f"GP{i}"] for i in range(60)})
SKEL = {"1": "SK_A", "2": "SK_B", "3": "SK_C", "4": "SK_H", "5": "SK_E", "6": "SK_G"}
feats = [
    {"label": "cmpA", "chebi": {"1"}, "skeleton": "SK_A"},                    # exact
    {"label": "cmpB", "chebi": {"99"}, "skeleton": "SK_B"},                    # not in file by ID; same skeleton as 2
    {"label": "cmpHub", "chebi": {"4"}, "skeleton": "SK_H"},                   # 61 reactions > 60
    {"label": "cmpEdge", "chebi": {"6"}, "skeleton": "SK_G"},                  # exactly 60 reactions: allowed
    {"label": "cmpNone", "chebi": {"77"}, "skeleton": "SK_Z"},                 # nothing in Reactome
    {"label": "cmpNoProt", "chebi": {"5"}, "skeleton": "SK_E"},                # reaction lists no protein
]
res = mn.metabolite_enzymes(feats, CHEBI_RX, RX_PROT, SKEL, 60)
st = {r["Compound"]: r for r in res["rows"]}
check(st["cmpA"]["Status"] == "used" and st["cmpA"]["Reactions"] == 2 and st["cmpA"]["Proteins"] == 2 and st["cmpA"]["Matched via"] == "ChEBI ID",
      "an exact ChEBI match seeds the proteins of its reactions")
check(st["cmpB"]["Status"] == "used" and "skeleton" in st["cmpB"]["Matched via"] and st["cmpB"]["Matched ChEBI"] == "2" and st["cmpB"]["Proteins"] == 3,
      "SENSITIVITY: a compound Reactome cites under another ChEBI number is matched on its InChIKey skeleton")
check(st["cmpHub"]["Status"].startswith("left out: in 61 reactions") and st["cmpEdge"]["Status"] == "used",
      "SENSITIVITY: the hub guard excludes 61 reactions and keeps exactly 60 (boundary)")
check(st["cmpNone"]["Status"].startswith("left out: not a participant") and st["cmpNoProt"]["Status"].startswith("left out: its reactions list no protein"),
      "a compound with no reaction, and one whose reactions have no protein, are left out WITH the reason")
check(set(res["proteins"]) == {"P1", "P2", "P3"} | {f"GP{i}" for i in range(60)} and res["proteins"]["P2"] == ["cmpA", "cmpB"],
      "the protein table says which compounds seeded each protein (P2 is in the reactions of both A and B)")
res_low = mn.metabolite_enzymes(feats, CHEBI_RX, RX_PROT, SKEL, 59)
check({r["Compound"]: r["Status"] for r in res_low["rows"]}["cmpEdge"].startswith("left out: in 60 reactions"), "lowering the limit to 59 excludes the 60-reaction compound")

# the second guard: a reaction that lists a whole family of proteins (a ligand-binding step lists every receptor)
FAM_RX = dict(CHEBI_RX)
FAM_RX.update({"7": ["F30", "F31", "S1"], "8": ["F31"], "9": ["F30"]})
FAM_PROT = dict(RX_PROT)
FAM_PROT.update({"F30": [f"FP{i}" for i in range(30)], "F31": [f"FQ{i}" for i in range(31)], "S1": ["PS1", "PS2"]})
FAM_SKEL = dict(SKEL)
fam = mn.metabolite_enzymes([{"label": "ligand", "chebi": {"7"}, "skeleton": ""}, {"label": "onlyFamily", "chebi": {"8"}, "skeleton": ""},
                             {"label": "exactly30", "chebi": {"9"}, "skeleton": ""}], FAM_RX, FAM_PROT, FAM_SKEL, 60, 30)
fs = {r["Compound"]: r for r in fam["rows"]}
check(fs["ligand"]["Status"] == "used" and fs["ligand"]["Family reactions skipped"] == 1 and fs["ligand"]["Proteins"] == 32 and "FQ0" not in fam["proteins"],
      "SENSITIVITY: a reaction listing 31 proteins is skipped (counted), the 30-protein and the 2-protein reactions of the same compound still seed")
check(fs["onlyFamily"]["Status"].startswith("left out: its reactions list no protein, or only large protein families"),
      "a compound whose ONLY reactions are family reactions is left out, and the reason says so")
check(fs["exactly30"]["Status"] == "used" and fs["exactly30"]["Proteins"] == 30, "BOUNDARY: a reaction with exactly 30 proteins is kept")
fam_hi = mn.metabolite_enzymes([{"label": "onlyFamily", "chebi": {"8"}, "skeleton": ""}], FAM_RX, FAM_PROT, FAM_SKEL, 60, 31)
check(fam_hi["rows"][0]["Status"] == "used" and fam_hi["rows"][0]["Proteins"] == 31, "raising the protein limit to 31 lets the 31-protein reaction seed")


def ref_enzymes(features, chebi_rx, rx_prot, skeletons, cap, pcap=10 ** 9):
    """Independent: loop over the FILE's ChEBI ids and test each against the compound (no skeleton index)."""
    out = {}
    for f in features:
        rx = set()
        for c in chebi_rx:
            if c in f["chebi"] or (f["skeleton"] and skeletons.get(c) == f["skeleton"]):
                rx |= set(chebi_rx[c])
        prot = set()
        for r in rx:
            if len(rx_prot.get(r, [])) <= pcap:
                prot |= set(rx_prot.get(r, []))
        out[f["label"]] = (len(rx), prot if (rx and len(rx) <= cap and prot) else None)
    return out


rg = np.random.RandomState(3)
chebi_big = {str(i): [f"R{rg.randint(0, 60)}" for _ in range(rg.randint(1, 8))] for i in range(80)}
rx_big = {f"R{i}": [f"P{rg.randint(0, 40)}" for _ in range(rg.randint(0, 4))] for i in range(60)}
skel_big = {str(i): f"S{rg.randint(0, 50)}" for i in range(80)}
feat_big = [{"label": f"c{i}", "chebi": {str(rg.randint(0, 100))}, "skeleton": f"S{rg.randint(0, 60)}"} for i in range(60)]
for cap, pcap in ((3, 10 ** 9), (8, 2), (20, 3), (500, 1), (500, 10 ** 9)):
    got = mn.metabolite_enzymes(feat_big, chebi_big, rx_big, skel_big, cap, pcap)
    ref = ref_enzymes(feat_big, chebi_big, rx_big, skel_big, cap, pcap)
    ok = True
    for r in got["rows"]:
        n_rx, prot = ref[r["Compound"]]
        used = r["Status"] == "used"
        ok &= (r["Reactions"] == n_rx) and (used == (prot is not None))
        if used:
            ok &= {a for a, labs in got["proteins"].items() if r["Compound"] in labs} == prot
    check(ok, f"60 random compounds against the brute-force join at limits of {cap} reactions / {pcap} proteins per reaction: statuses, counts and proteins agree")

# =============================================================================
print("=== Genomics seeds: UniProt first, symbol second, the rest counted ===")
acc_idx = {"by_accession": {"ACC1": "p1", "ACC4": "p4"}, "by_symbol": {"GENE2": "p2", "GENE1": "p1x"}}
gdf = pd.DataFrame({"UniProt ID": ["ACC1", "", "NOPE", "ACC4", ""], "Gene Symbol": ["GENE1", "GENE2", "GENE9", "GENE4", ""]})
gs = mn.genomics_seeds(gdf, acc_idx)
check(gs["seeds"] == {"p1", "p2", "p4"} and gs["n_unmapped"] == 2 and gs["n_input"] == 5,
      f"SENSITIVITY: the UniProt ID wins (ACC1 -> p1, not GENE1's p1x), the symbol is the fallback (GENE2 -> p2), the rest are counted ({sorted(gs['seeds'])}, {gs['n_unmapped']} unmapped)")

print("=== a STRING-like fixture: 3 communities of 20 proteins ===")
# =============================================================================
TAX = "9995"
N = 60
names = [f"Q{i}" for i in range(N)]
PID = {n: f"{TAX}.{n}" for n in names}
rng = np.random.default_rng(5)
edges = []
for c in range(3):
    base = c * 20
    for i in range(20):
        for j in range(i + 1, 20):
            if rng.random() < 0.45:
                edges.append((names[base + i], names[base + j], int(rng.integers(750, 999))))
    for i in range(19):
        edges.append((names[base + i], names[base + i + 1], 900))
for a, b in [(5, 25), (30, 45), (10, 50)]:
    edges.append((names[a], names[b], 720))
links = []
for a, b, w in edges:
    links += [(PID[a], PID[b], w), (PID[b], PID[a], w)]


def _gz(path, header, rows, sep):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        fh.write(sep.join(header) + "\n")
        for r in rows:
            fh.write(sep.join(str(x) for x in r) + "\n")


TMP = tempfile.TemporaryDirectory()
LINKS, INFO, ALIASES = (os.path.join(TMP.name, f) for f in ("l.gz", "i.gz", "a.gz"))
_gz(LINKS, ["protein1", "protein2", "combined_score"], links, " ")
_gz(INFO, ["#string_protein_id", "preferred_name", "protein_size", "annotation"], [(PID[n], f"GENE{i}", "100", "x") for i, n in enumerate(names)], "\t")
_gz(ALIASES, ["#string_protein_id", "alias", "source"], [(PID[n], f"ACC{i}", "UniProt_AC") for i, n in enumerate(names)], "\t")
_real_paths, _real_in_string = gnet._string_file_paths, gnet.is_taxon_in_string
gnet._string_file_paths = lambda taxon_id: {"links": LINKS, "info": INFO, "aliases": ALIASES}
gnet.is_taxon_in_string = lambda taxon_id: True
FAKE_GO = {"biological_process": {"GO:0001": {"name": "community A process", "genes": {f"ACC{i}" for i in range(0, 20)}},
                                   "GO:0002": {"name": "community B process", "genes": {f"ACC{i}" for i in range(20, 40)}}},
           "molecular_function": {}, "cellular_component": {}}
_real_go = gm.build_go_universe
gm.build_go_universe = lambda *a, **k: FAKE_GO


def prot_df(sig):
    rows = []
    for i in ROW_ORDER:
        s = i in sig
        rows.append({"input_id": f"r{i}", "accession": f"ACC{i}", "gene_symbol": f"GENE{i}", "protein_name": "",
                     "log2FC": 3.0 if s else 0.05, "pvalue": 1e-8 if s else 0.8, "padj": 1e-6 if s else 0.9, "status": "resolved"})
    return pd.DataFrame(rows)


def tx_df(sig):
    rows = []
    for i in ROW_ORDER:
        s = i in sig
        rows.append({"input_id": f"t{i}", "gene_symbol": f"GENE{i}", "ensembl_id": "", "entrez_id": "",
                     "log2FC": -2.5 if s else 0.1, "pvalue": 1e-7 if s else 0.7, "padj": 1e-5 if s else 0.8, "status": "resolved"})
    return pd.DataFrame(rows)


# a 60-protein graph and 30 permutations cannot reach FDR < 0.05; a loose alpha gives hits (and so modules) to test the plumbing
CFG = {"n_permutations": 30, "confidence": 700, "min_module_size": 3, "seed": 42, "fdr_alpha": 0.5}
# The permutation null draws its random numbers seed by seed, so the existing tile (seeds in the order of the DE table's rows) and
# this engine (seeds in graph-index order) only coincide when the table's rows are in graph-index order; the fixture is built so.
_G = gnet.build_string_network(TAX, confidence=700)
ROW_ORDER = sorted(range(N), key=lambda i: _G["protein_index"][PID[f"Q{i}"]])

try:
    # =========================================================================
    print("=== SAME SEEDS -> the existing Proteomics tile, exactly ===")
    # =========================================================================
    sig = set(range(0, 6)) | {22, 25, 41, 44}
    pdf = prot_df(sig)
    single = pnet.run_ppi_network_analysis(pdf, TAX, taxon_name="Test organism", confidence=700, n_permutations=30, min_module_size=3, seed=42, fdr_alpha=0.5)
    # Genomics lists the very same proteins: both omics have identical seed sets, so every weight is exactly 1.0
    gen_df = pd.DataFrame({"UniProt ID": [f"ACC{i}" for i in sorted(sig)], "Gene Symbol": [f"GENE{i}" for i in sorted(sig)]})
    both = mn.run_string_network({"genomics": {"final_df": gen_df}, "proteomics": {"resolved_df": pdf}}, TAX, "Test organism", CFG)
    a = single["propagation"]["table"].set_index("STRING ID")
    b = both["propagation"]["table"].set_index("STRING ID")
    cols = ["Propagated score", "Z-score", "p-value", "Empirical p-value", "FDR", "FDR (recruited family)", "Seed score"]
    diffs = {c: float((a.loc[b.index, c].astype(float).fillna(-1) - b[c].astype(float).fillna(-1)).abs().max()) for c in cols}
    same = all(d < 1e-9 for d in diffs.values())
    check(same and set(a.index) == set(b.index) and (a.loc[b.index, "Seed"] == b["Seed"]).all(),
          f"THE KEY CHECK: identical seeds from two omics give the existing Proteomics network tile's scores, Z, p, FDR and seeds (largest differences {diffs})")
    rwr_same = np.allclose(a.loc[b.index, "Propagated score"].astype(float), b["Propagated score"].astype(float), rtol=1e-12)
    check(rwr_same, "the random-walk scores themselves (no random numbers involved) equal the tile's")
    check(both["propagation"]["n_hits"] == single["propagation"]["n_hits"] and both["propagation"]["n_recruited"] == single["propagation"]["n_recruited"],
          f"...and the same number of hits ({both['propagation']['n_hits']}) and recruited proteins ({both['propagation']['n_recruited']})")
    check(both["modules"] is not None and single["modules"] is not None and
          sorted(map(sorted, both["modules"]["membership"].values())) == sorted(map(sorted, single["modules"]["membership"].values())),
          f"...and the same modules ({0 if both['modules'] is None else len(both['modules']['membership'])})")

    # =========================================================================
    print("=== different seeds from three omics ===")
    # =========================================================================
    p_sig, t_sig = set(range(0, 6)), set(range(3, 9)) | {22}
    enzyme_maps = {"status": "ok", "chebi_reactions": {"10": ["RA"], "11": ["RB"]}, "reaction_proteins": {"RA": ["ACC40", "ACC41", "ACC42"], "RB": ["ACC43"]}}
    feats3 = [{"label": "metabA", "chebi": {"10"}, "skeleton": ""}, {"label": "metabB", "chebi": {"11"}, "skeleton": ""}]
    res3 = mn.run_string_network(
        {"proteomics": {"resolved_df": prot_df(p_sig)}, "transcriptomics": {"resolved_df": tx_df(t_sig)}, "metabolomics": {"features": feats3}},
        TAX, "Test organism", CFG, enzyme_maps=enzyme_maps, skeleton_provider=lambda ids: {})
    tbl = res3["propagation"]["table"].set_index("Gene Symbol")
    check(res3["omics_used"] == ["transcriptomics", "proteomics", "metabolomics"], f"three omics used, in the fixed order ({res3['omics_used']})")
    check(tbl.loc["GENE4", "Seed omics"] == "Transcriptomics + Proteomics" and tbl.loc["GENE0", "Seed omics"] == "Proteomics"
          and tbl.loc["GENE8", "Seed omics"] == "Transcriptomics" and tbl.loc["GENE41", "Seed omics"] == "Metabolomics",
          "each seed names the omics that seeded it (GENE4 is hit by transcripts AND proteins; GENE41 only through a metabolite's reaction)")
    check(tbl.loc["GENE4", "Seed score"] > tbl.loc["GENE0", "Seed score"] > 0, "a protein hit by two omics has the larger seed score")
    n_seeds = res3["seeds"]["n_seeds"]
    check(n_seeds == len(p_sig | t_sig) + 4 and abs(tbl["Seed score"].dropna().sum() - n_seeds) < 1e-9, f"{n_seeds} distinct seeds; their weights sum to that number")
    check(abs(sum(r["Share of the combined weight"] for r in res3["omics"]) - 1.0) < 1e-12, "the per-omics shares of the combined weight sum to 1")
    check(abs(float(tbl.loc["GENE40", "Transcriptomics log2FC"]) - 0.1) < 1e-12 and tbl.loc["GENE40", "Seed omics"] == "Metabolomics",
          "a metabolite-seeded protein still shows the transcript level the user measured for it (not significant: log2FC 0.1)")
    check(abs(float(tbl.loc["GENE4", "Transcriptomics log2FC"]) + 2.5) < 1e-12 and abs(float(tbl.loc["GENE4", "Proteomics log2FC"]) - 3.0) < 1e-12,
          "the user's own measurements travel with the protein, each under its own omics")
    check({"Seeds from Transcriptomics", "Seeds from Proteomics", "Seeds from Metabolomics", "Seeds from several omics"} <= set(res3["modules"]["table"].columns),
          "the module table counts the seeds each omics contributed")
    mt = res3["modules"]["table"]
    check((mt["Seeds"] == mt[["Seeds from Transcriptomics", "Seeds from Proteomics", "Seeds from Metabolomics"]].sum(axis=1) - mt["Seeds from several omics"]).all(),
          "per module: seeds = sum of per-omics seeds minus the ones hit by several omics (no double count)")
    check(res3["metabolites"] is not None and {r["Compound"] for r in res3["metabolites"]} == {"metabA", "metabB"}, "the metabolite -> enzyme table is returned")
    check(any("association" in c for c in res3["caveats"]) and any("Equal weight per omics" in c for c in res3["caveats"]),
          "the caveats state the equal-weight rule and that metabolite seeds are an association")
    check(res3["propagation"]["n_proteins"] == N and "_graph" in res3, "the payload carries the graph in-process like the other network tiles")
    # the permutation null draws its random numbers seed by seed, so the ORDER seeds enter it matters; a Python set has no stable
    # order, so the engine must hand them over in graph-index order
    _orig_core = gnet._propagation_core
    _seen = {}

    def _spy(graph, seed_idx_score, *a, **k):
        _seen["order"] = list(seed_idx_score)
        return _orig_core(graph, seed_idx_score, *a, **k)
    gnet._propagation_core = _spy
    try:
        mn.run_string_network({"proteomics": {"resolved_df": prot_df(p_sig)}, "transcriptomics": {"resolved_df": tx_df(t_sig)}}, TAX, "Test organism", CFG)
    finally:
        gnet._propagation_core = _orig_core
    check(len(_seen["order"]) > 8 and _seen["order"] == sorted(_seen["order"]),
          f"SENSITIVITY: the seeds enter the permutation null in graph-index order, not in a set's arbitrary order ({_seen['order'][:6]}...)")

    # =========================================================================
    print("=== omics that cannot take part ===")
    # =========================================================================
    name, msg = raises(mn.run_string_network, {"proteomics": {"resolved_df": prot_df(p_sig)}}, TAX, "Test organism", CFG)
    check(name == "ValueError" and "at least two omics" in msg, f"one omics alone is refused ({msg[:60]})")
    none_sig = prot_df(set())
    name, msg = raises(mn.run_string_network, {"proteomics": {"resolved_df": none_sig}, "transcriptomics": {"resolved_df": tx_df(t_sig)}}, TAX, "Test organism", CFG)
    check(name == "ValueError" and "at least two omics" in msg and "Proteomics has no seed" in msg, "an omics whose DE table has no significant row is left out, with the reason")
    name, msg = raises(mn.run_string_network, {"proteomics": {"resolved_df": prot_df(p_sig)}, "metabolomics": {"features": feats3}}, TAX, "Test organism", CFG,
                       enzyme_maps={"status": "error", "note": "no network"})
    check(name == "ValueError" and "Metabolomics was left out" in msg and "no network" in msg, "unavailable Reactome files leave metabolomics out and say why")
finally:
    gnet._string_file_paths, gnet.is_taxon_in_string = _real_paths, _real_in_string
    gm.build_go_universe = _real_go

# =============================================================================
print("=== the ChEBI skeleton cache (negatives remembered, failures not) ===")
import Multiomics_module_retriever as mr  # noqa: E402

_store = {}
_orig = (mr._gm._load_universe_cache, mr._gm._save_universe_cache, mr.fetch_chebi_skeletons)
mr._gm._load_universe_cache = lambda db, key: _store.get((db, key))
mr._gm._save_universe_cache = lambda db, key, data: _store.__setitem__((db, key), data)
_calls = []


def _fake_fetch(answer):
    def f(ids):
        _calls.append(sorted(ids))
        return {i: "SK" + i for i in ids if i in answer}
    return f


try:
    ids10 = {str(i) for i in range(10)}
    mr.fetch_chebi_skeletons = _fake_fetch({str(i) for i in range(9)})        # ChEBI has no structure for id 9
    got = mr.fetch_chebi_skeletons_cached(ids10)
    check(set(got) == {str(i) for i in range(9)} and got["3"] == "SK3" and len(_calls) == 1, "first call asks ChEBI once and returns the ids that have a structure")
    got2 = mr.fetch_chebi_skeletons_cached(ids10)
    check(got2 == got and len(_calls) == 1, "SENSITIVITY: the second call asks nobody - the id with no structure is remembered, so it is not asked for again")
    got3 = mr.fetch_chebi_skeletons_cached(ids10 | {"20", "21"})
    check(_calls[-1] == ["20", "21"] and len(_calls) == 2 and set(got3) == set(got), "a new call asks only for the ids it has never seen")
    _store.clear()
    _calls.clear()
    mr.fetch_chebi_skeletons = _fake_fetch({"0", "1"})                         # only 2 of 10 answered: the lookup mostly failed
    got4 = mr.fetch_chebi_skeletons_cached(ids10)
    mr.fetch_chebi_skeletons_cached(ids10)
    check(set(got4) == {"0", "1"} and len(_calls) == 2 and _calls[1] == sorted(ids10 - {"0", "1"}),
          "SENSITIVITY: when most ids did not answer (a failed batch) the silence is NOT remembered as 'no structure' - the next call asks again")
    _store.clear()
    _calls.clear()
    mr.fetch_chebi_skeletons = _fake_fetch(set())
    check(mr.fetch_chebi_skeletons_cached(ids10) == {} and not _store and mr.fetch_chebi_skeletons_cached(set()) == {}, "a total failure returns nothing and stores nothing")
finally:
    mr._gm._load_universe_cache, mr._gm._save_universe_cache, mr.fetch_chebi_skeletons = _orig

print("=== Reactome's reaction files: a flaky download is retried, a dead one is reported ===")
_orig_stream, _orig_sleep = mr._stream_reactome_pairs, mr._time.sleep
mr._time.sleep = lambda s: None
_attempts = {"n": 0}


def _flaky(url, species):
    _attempts["n"] += 1
    if _attempts["n"] <= 2:
        yield ("1", "R1", "a reaction")              # half the file arrives ...
        raise ConnectionError("Response ended prematurely")   # ... and the connection drops
    if "ChEBI" in url:
        yield ("16947", "R1", "a reaction")
    else:
        yield ("P12345", "R1", "a reaction")


mr._stream_reactome_pairs = _flaky
try:
    rows_ok = mr._read_reactome_rows("https://x/ChEBI2ReactomeReactions.txt", "homo sapiens")
    check(rows_ok == [("16947", "R1", "a reaction")] and _attempts["n"] == 3,
          "SENSITIVITY: a download that drops twice is read again from the start, and the rows of the failed attempts are NOT kept")


    def _dead(url, species):
        raise ConnectionError("HTTP 522")
        yield
    mr._stream_reactome_pairs = _dead
    name_e, msg_e = raises(mr._read_reactome_rows, "https://x/file.txt", "homo sapiens")
    check(name_e == "ConnectionError" and "522" in msg_e, "a download that never works raises after the attempts")
    _cache_store = {}
    _o = (mr._gm._load_universe_cache, mr._gm._save_universe_cache)
    mr._gm._load_universe_cache = lambda db, key: _cache_store.get((db, key))
    mr._gm._save_universe_cache = lambda db, key, data: _cache_store.__setitem__((db, key), data)
    try:
        res_dead = mr.fetch_reactome_reaction_maps("9606", "Homo sapiens")
    finally:
        mr._gm._load_universe_cache, mr._gm._save_universe_cache = _o
    check(res_dead["status"] == "error" and "Try again" in res_dead["note"] and not _cache_store, "the app reports an error (and caches nothing) when Reactome stays down")
finally:
    mr._stream_reactome_pairs, mr._time.sleep = _orig_stream, _orig_sleep

print("=== the Reactome pathway-overlap network ===")
# =============================================================================
TERMS = {
    "T1": {"name": "Alpha", "genes": {"a1", "a2", "a3", "a4"}},
    "T2": {"name": "Beta", "genes": {"a3", "a4", "a5", "a6"}},          # shares a3,a4 with T1 -> hits {A3,A4} vs ...
    "T3": {"name": "Gamma", "genes": {"z1", "z2"}},
    "T4": {"name": "Delta (metabolites)", "genes": {"m1", "m2"}},
    "T5": {"name": "Epsilon, no hits", "genes": {"q1"}},
}
inputs = {
    "transcriptomics": {"term_map": TERMS, "list_set": {"a1", "a2", "a3", "a4", "a5", "z1"}, "universe_set": set(), "labels": None},
    "metabolomics": {"term_map": {"T4": TERMS["T4"], "T1": {"name": "Alpha", "genes": {"m1"}}}, "list_set": {"m1", "m2"}, "universe_set": set(),
                     "labels": {"m1": "A3", "m2": "Zeta"}},       # a compound NAMED exactly like gene A3 must not merge with it
}
block = {"table": pd.DataFrame({"Term ID": ["T1", "T2", "T3", "T4", "T5"], "Term Name": ["Alpha", "Beta", "Gamma", "Delta", "Epsilon"],
                                "Combined p-value": [0.001, 0.002, 0.01, 0.02, 0.5], "Combined FDR": [0.01, 0.02, 0.04, 0.2, 0.9]}),
         "omics_combined": ["transcriptomics", "metabolomics"], "note": None}
net = mn.pathway_overlap_network(block, inputs, top_n=10, min_jaccard=0.25)
ids = [n["id"] for n in net["nodes"]]
check(ids == ["T1", "T2", "T3", "T4"], f"nodes are the pathways with a hit, best combined p first - T5 (no hit) is dropped ({ids})")
pool = {"T1": {"A1", "A2", "A3", "A4", "M:A3"}, "T2": {"A3", "A4", "A5"}, "T3": {"Z1"}, "T4": {"M:Zeta"}}
ref_edges = {}
for i in pool:
    for j in pool:
        if i < j and pool[i] & pool[j]:
            ref_edges[(i, j)] = len(pool[i] & pool[j]) / len(pool[i] | pool[j])
got_edges = {(e["a"], e["b"]): e["jaccard"] for e in net["edges"]}
ref_keep = {k: v for k, v in ref_edges.items() if v >= 0.25}
check(got_edges.keys() == ref_keep.keys() and all(abs(got_edges[k] - ref_keep[k]) < 1e-12 for k in ref_keep),
      f"edges equal the brute-force Jaccard of the pooled hit members ({ {k: round(v, 3) for k, v in got_edges.items()} })")
check(abs(ref_edges[("T1", "T2")] - 2 / 6) < 1e-12, "(T1, T2) share A3 and A4 of 6 distinct members = 1/3, above the 0.25 cut")
net_hi = mn.pathway_overlap_network(block, inputs, top_n=10, min_jaccard=1 / 3)
net_hi2 = mn.pathway_overlap_network(block, inputs, top_n=10, min_jaccard=1 / 3 + 1e-9)
check(len(net_hi["edges"]) == 1 and len(net_hi2["edges"]) == 0, "SENSITIVITY: the threshold boundary is inclusive (an edge at exactly 1/3 stays at 1/3, goes at 1/3 + 1e-9)")
check(net["nodes"][0]["n_hit_members"] == 5,
      "SENSITIVITY: a compound named exactly like a gene (A3) is NOT merged with the gene (T1 has 5 distinct hit members, not 4)")
grp = {g["id"]: set(g["terms"]) for g in net["groups"]}
check(grp[1] == {"T1", "T2"} and {frozenset(v) for v in grp.values()} == {frozenset({"T1", "T2"}), frozenset({"T3"}), frozenset({"T4"})},
      "groups are the connected components (T1-T2 together; T3 and T4 alone)")
check(net["nodes"][0]["significant"] and not net["nodes"][3]["significant"] and net["nodes"][0]["hits"] == {"transcriptomics": 4, "metabolomics": 1},
      "nodes carry the combined FDR flag and the hits per omics")
check(len(mn.pathway_overlap_network(block, inputs, top_n=2)["nodes"]) == 2, "top_n limits the nodes (by combined p-value)")
empty = mn.pathway_overlap_network({"table": pd.DataFrame(), "omics_combined": [], "note": "x"}, inputs)
check(empty["nodes"] == [] and "no pathway network" in empty["note"].lower(), "an empty joint result gives a note, not an error")

print("=== near-duplicate pathways are merged into the better-ranked one ===")
TERMS2 = dict(TERMS)
TERMS2["T6"] = {"name": "Alpha sub-step", "genes": {"a1", "a2", "a3", "a4"}}        # hit members {A1..A4}; T1's pool is {A1..A4, M:A3}: Jaccard 4/5
inputs2 = {"transcriptomics": dict(inputs["transcriptomics"], term_map=TERMS2), "metabolomics": inputs["metabolomics"]}
block2 = {"table": pd.DataFrame({"Term ID": ["T1", "T6", "T2", "T3", "T4"], "Term Name": ["Alpha", "Alpha sub-step", "Beta", "Gamma", "Delta"],
                                 "Combined p-value": [0.001, 0.0015, 0.002, 0.01, 0.02], "Combined FDR": [0.01, 0.012, 0.02, 0.04, 0.2]}),
          "omics_combined": ["transcriptomics", "metabolomics"], "note": None}
m1 = mn.pathway_overlap_network(block2, inputs2, top_n=10, min_jaccard=0.25, merge_jaccard=0.7)
check([n["id"] for n in m1["nodes"]] == ["T1", "T2", "T3", "T4"] and m1["nodes"][0]["merged"] == [{"id": "T6", "name": "Alpha sub-step"}] and m1["n_merged"] == 1,
      "SENSITIVITY: T6 (Jaccard 0.8 with T1) is merged into T1, the better-ranked one, and listed on it")
m_edge = mn.pathway_overlap_network(block2, inputs2, top_n=10, merge_jaccard=0.8)
m_above = mn.pathway_overlap_network(block2, inputs2, top_n=10, merge_jaccard=0.8 + 1e-9)
check(m_edge["n_merged"] == 1 and m_above["n_merged"] == 0 and "T6" in [n["id"] for n in m_above["nodes"]], "the merge threshold is inclusive (merged at exactly 0.8, kept apart at 0.8 + 1e-9)")
m_n = mn.pathway_overlap_network(block2, inputs2, top_n=3, merge_jaccard=0.7)
check([n["id"] for n in m_n["nodes"]] == ["T1", "T2", "T3"], "SENSITIVITY: top_n counts DISTINCT pathways - the merged T6 does not use up a place")
check(not any(abs(len(a & b) / len(a | b)) >= 0.7 for a, b in [({"A1"}, {"A9"})]), "(sanity)")
ids_drawn = [n["id"] for n in m1["nodes"]]
pool2 = {"T1": {"A1", "A2", "A3", "A4", "M:A3"}, "T2": {"A3", "A4", "A5"}, "T3": {"Z1"}, "T4": {"M:Zeta"}}
check(all(len(pool2[a] & pool2[b]) / len(pool2[a] | pool2[b]) < 0.7 for i, a in enumerate(ids_drawn) for b in ids_drawn[i + 1:]), "no two drawn pathways overlap at or above the merge threshold")
check(all(e["jaccard"] < 0.7 for e in m1["edges"]), "so every link is between the link threshold and the merge threshold")

print("=== layout ===")
lid = [f"n{i}" for i in range(40)]
lr = np.random.RandomState(1)
ledges = [(lid[i], lid[j], float(lr.rand())) for i in range(40) for j in range(i + 1, 40) if lr.rand() < 0.08]
pos1, pos2 = mn.force_layout(lid, ledges), mn.force_layout(lid, ledges)
P1 = np.array([pos1[t] for t in lid])
dm = np.sqrt(((P1[:, None, :] - P1[None, :, :]) ** 2).sum(axis=2)) + np.eye(40) * 9
check(pos1 == pos2, "the layout is deterministic")
check(dm.min() >= 0.09 - 1e-6, f"SENSITIVITY: no two of 40 nodes are closer than the minimum gap (closest {dm.min():.3f})")
adj = np.zeros((40, 40), bool)
for a, b, _ in ledges:
    adj[lid.index(a), lid.index(b)] = adj[lid.index(b), lid.index(a)] = True
check(dm[adj].mean() < dm[~adj & ~np.eye(40, dtype=bool)].mean(), "connected nodes sit closer together than unconnected ones on average")
check(mn.force_layout([], []) == {} and mn.force_layout(["x"], []) == {"x": (0.5, 0.5)}, "0 and 1 nodes are handled")

print("=== the pathway canvas: label boxes never overlap ===")
BW, BH = mn.LABEL_BOX
for n_nodes, dens in ((12, 0.2), (35, 0.08), (60, 0.05), (100, 0.03)):
    ids_c = [f"p{i}" for i in range(n_nodes)]
    lr2 = np.random.RandomState(n_nodes)
    ed_c = [(ids_c[i], ids_c[j], float(lr2.rand())) for i in range(n_nodes) for j in range(i + 1, n_nodes) if lr2.rand() < dens]
    posc, wc, hc = mn.pathway_canvas_layout(ids_c, ed_c)
    pts = np.array([posc[t] for t in ids_c])
    worst = max((min(BW - abs(pts[i, 0] - pts[j, 0]), BH - abs(pts[i, 1] - pts[j, 1])) for i in range(n_nodes) for j in range(i + 1, n_nodes)), default=-1)
    inside = (pts[:, 0] - BW / 2 >= -1e-6).all() and (pts[:, 1] - BH / 2 >= -1e-6).all() and (pts[:, 0] + BW / 2 <= wc + 1e-6).all() and (pts[:, 1] + BH / 2 <= hc + 1e-6).all()
    check(worst <= 1e-6 and inside, f"SENSITIVITY: {n_nodes} nodes - no two {BW:.0f}x{BH:.0f} label boxes overlap (brute force over {n_nodes * (n_nodes - 1) // 2} pairs) and every box is inside the {wc}x{hc} canvas")
a_pos, a_w, a_h = mn.pathway_canvas_layout(ids_c, ed_c)
b_pos, b_w, b_h = mn.pathway_canvas_layout(ids_c, ed_c)
check(a_pos == b_pos and (a_w, a_h) == (b_w, b_h), "the canvas layout is deterministic")
check(mn.pathway_canvas_layout([], []) == ({}, 1100.0, 800.0) and len(mn.pathway_canvas_layout(["a"], [])[0]) == 1, "0 and 1 pathways are handled")
two = np.array([[100.0, 100.0], [110.0, 400.0]])
mn.separate_boxes(two, BW, BH)
check(abs(two[0, 0] - two[1, 0]) < 1e-9 + 10 and abs(two[0, 1] - 100) < 1e-9, "boxes that already clear each other vertically are left alone")
close = np.array([[100.0, 100.0], [180.0, 110.0]])
mn.separate_boxes(close, BW, BH)
check(abs(close[0, 0] - close[1, 0]) >= BW - 1e-6 or abs(close[0, 1] - close[1, 1]) >= BH - 1e-6, "two overlapping boxes are pushed apart along the shorter way out")

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
