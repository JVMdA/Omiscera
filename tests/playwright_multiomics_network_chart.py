"""Fast browser test of the Network-based Integration results screen on a SYNTHETIC payload (seconds, no run, no network).

The real run on the bundled examples finds no module and no recruited protein (it is a small effect on a diffuse seed set), so
this test builds a payload that HAS them - modules with seeds from several omics, recruited proteins, per-module enrichment, a
metabolite table and a crowded pathway network (~30 pathways with long Reactome-style names, laid out by the engine's own
force_layout) - and injects it into the page. It checks what only a drawn page can show:
  - the panes, the recruited-only default and its toggle, the per-omics module columns;
  - the module composition chart (one series per omics + recruited, stacked, bars add up to the module's own counts);
  - the single-omics tiles' figures (network, module network, null distribution, enrichment dot plot) draw from this payload;
  - the pathway network: every node drawn, the SAME node positions as the engine's, no two labels overlapping (measured), chips
    present only where an omics hit the pathway, significant pathways outlined, group colours, the tools row;
  - no page error.
Mutation-testable: each rule below goes red when the matching line of app.js is broken.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from playwright.sync_api import sync_playwright

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
import Multiomics_network as mn  # noqa: E402

BASE = "http://127.0.0.1:8009"
PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


# ------------------------------------------------------------------ the synthetic payload
rng = np.random.RandomState(4)
OM = ["genomics", "transcriptomics", "proteomics", "metabolomics"]
LAB = {"genomics": "Genomics", "transcriptomics": "Transcriptomics", "proteomics": "Proteomics", "metabolomics": "Metabolomics"}
NAMES = ["Ribosome-associated quality control", "Metabolism of amino acids and derivatives", "Signaling by ROBO receptors", "Translation",
         "Peptide chain elongation", "PELO:HBS1L and ABCE1 dissociate a ribosome", "Eukaryotic Translation Termination", "Selenoamino acid metabolism",
         "Nonsense-Mediated Decay (NMD)", "Citric acid cycle (TCA cycle)", "Pyruvate metabolism and Citric Acid (TCA) cycle", "Glycolysis",
         "Gluconeogenesis", "Metabolism of carbohydrates", "Respiratory electron transport", "Complex I biogenesis", "Fatty acid metabolism",
         "Mitochondrial protein import", "Cellular response to starvation", "Regulation of expression of SLITs and ROBOs", "Glucose metabolism",
         "Transport of small molecules", "Metabolism of lipids", "Branched-chain amino acid catabolism", "Amino acid transport across the plasma membrane",
         "SLC-mediated transmembrane transport", "Pentose phosphate pathway", "Urea cycle", "Metabolism of nucleotides", "Purine catabolism",
         "Biological oxidations", "Phase II - Conjugation of compounds", "Glutathione conjugation", "Aerobic respiration and respiratory electron transport",
         "Mitochondrial translation"]
genes = [f"G{i}" for i in range(60)]
terms, hit_t, hit_p, hit_g = {}, set(), set(), set()
base_sets = [set(rng.choice(60, size=rng.randint(6, 14), replace=False)) for _ in NAMES]
for k in range(1, 12):                                        # make clusters of overlapping pathways
    base_sets[k] = base_sets[k] | base_sets[0] if k % 3 == 0 else base_sets[k]
for i, nm in enumerate(NAMES):
    terms[f"R-HSA-{1000 + i}"] = {"name": nm, "genes": {genes[j] for j in base_sets[i]}}
for k in range(6):                                            # sub-steps with exactly the parent's members: Reactome has many such near-duplicates
    terms[f"R-HSA-{2000 + k}"] = {"name": NAMES[k] + " (sub-step)", "genes": set(terms[f"R-HSA-{1000 + k}"]["genes"])}
inp = {"transcriptomics": {"term_map": terms, "list_set": {genes[j] for j in range(0, 30)}, "universe_set": set(genes), "labels": None},
       "proteomics": {"term_map": terms, "list_set": {genes[j] for j in range(12, 26)}, "universe_set": set(genes), "labels": None},
       "genomics": {"term_map": terms, "list_set": {genes[j] for j in range(18, 24)}, "universe_set": set(genes), "labels": None}}
rows = []
for i, (tid, info) in enumerate(terms.items()):
    k_base = int(tid.split("-")[-1]) - (2000 if int(tid.split("-")[-1]) >= 2000 else 1000)
    pv = 10 ** -(30 - 0.8 * k_base) * (1.5 if tid.split("-")[-1].startswith("2") else 1.0)       # a sub-step ranks just after its parent
    rows.append({"Term ID": tid, "Term Name": info["name"], "Combined p-value": pv, "Combined FDR": min(1.0, pv * 40)})
block = {"table": pd.DataFrame(rows), "omics_combined": ["genomics", "transcriptomics", "proteomics"], "note": None}
for _thr in (0.25, 0.3, 0.35, 0.4, 0.45, 0.5):        # the lowest threshold that gives several groups and a real set of links
    pw = mn.pathway_overlap_network(block, inp, top_n=35, min_jaccard=_thr)
    if len(pw["groups"]) >= 4 and len(pw["edges"]) >= 8:
        break

# hit table: 40 proteins, 12 seeds from different omics, 8 recruited significant ones, 20 other
cols = mn.HIT_COLUMNS
hit_rows = []
seed_omics_of = {}
for i in range(40):
    is_seed = i < 14
    om = [["Genomics"], ["Transcriptomics"], ["Proteomics"], ["Transcriptomics", "Proteomics"], ["Metabolomics"]][i % 5] if is_seed else []
    sig = i < 22
    hit_rows.append({"STRING ID": f"9606.P{i}", "Gene Symbol": f"GENE{i}", "Mapped via": "", "log2FC": None, "padj": None, "Propagated score": 1e-3 / (i + 1),
                     "Z-score": 30 - i, "p-value": 1e-9 * (i + 1), "Empirical p-value": 0.002, "FDR": 1e-6 * (i + 1) if sig else 0.4,
                     "FDR (recruited family)": None if is_seed else (1e-5 * (i + 1) if sig else 0.5), "Seed": bool(is_seed),
                     "Seed score": 1.2 if is_seed else None, "Seed omics": " + ".join(om), "Transcriptomics log2FC": 2.0 if "Transcriptomics" in om else None,
                     "Transcriptomics padj": 1e-4 if "Transcriptomics" in om else None, "Proteomics log2FC": 1.5 if "Proteomics" in om else None, "Proteomics padj": 1e-3 if "Proteomics" in om else None})
edges = [{"source": f"9606.P{i}", "target": f"9606.P{j}", "weight": 0.9} for i in range(22) for j in range(i + 1, 22) if (i * 7 + j * 3) % 5 == 0]
mod_rows = [
    {"Module ID": "M1", "Size": 9, "Seeds": 6, "Recruited": 3, "Seeds from Genomics": 2, "Seeds from Transcriptomics": 4, "Seeds from Proteomics": 3, "Seeds from Metabolomics": 0,
     "Seeds from several omics": 3, "Genes": "A, B", "Top enriched term": "TCA cycle", "Top enriched FDR": 1e-5, "Top enriched FDR (raw, uncorrected)": 1e-6, "Enrichment database": "biological_process"},
    {"Module ID": "M2", "Size": 6, "Seeds": 4, "Recruited": 2, "Seeds from Genomics": 0, "Seeds from Transcriptomics": 1, "Seeds from Proteomics": 0, "Seeds from Metabolomics": 4,
     "Seeds from several omics": 1, "Genes": "C, D", "Top enriched term": "Urea cycle", "Top enriched FDR": 3e-3, "Top enriched FDR (raw, uncorrected)": 1e-3, "Enrichment database": "biological_process"},
]
enr_cols = ["Term ID", "Term Name", "k", "K", "p-value", "FDR", "Bonferroni", "Fold Enrichment", "Overlapping Genes"]
enr_tbl = {"columns": enr_cols, "rows": [{"Term ID": f"GO:{i}", "Term Name": f"term {i}", "k": 5 - i % 3, "K": 40, "p-value": 1e-4 * (i + 1), "FDR": 1e-3 * (i + 1), "Bonferroni": 0.01,
                                          "Fold Enrichment": 12.0 - i, "Overlapping Genes": "A, B"} for i in range(6)]}
metab = {"columns": ["Compound", "Matched ChEBI", "Matched via", "Reactions", "Family reactions skipped", "Proteins", "Status"],
         "rows": [{"Compound": "Citrate", "Matched ChEBI": "16947", "Matched via": "ChEBI ID", "Reactions": 12, "Family reactions skipped": 0, "Proteins": 15, "Status": "used"},
                  {"Compound": "ATP", "Matched ChEBI": "30616", "Matched via": "ChEBI ID", "Reactions": 2103, "Family reactions skipped": 300, "Proteins": 900,
                   "Status": "left out: in 2103 reactions (more than the limit of 60) - a hub such as a cofactor"}]}
payload = {
    "taxon_id": "9606", "taxon_name": "Homo sapiens",
    "parameters": {"string_version": "12.0", "confidence": 700, "restart_prob": 0.5, "n_permutations": 500, "permutation_type": "degree-preserving", "seed_rng": 42, "fdr_alpha": 0.05,
                   "mcl_inflation": 2.0, "min_module_size": 3, "max_hits_for_modules": 600, "enrichment_databases": ["GO"], "max_reactions": 60, "max_reaction_proteins": 30,
                   "weighting": "equal weight per omics"},
    "omics": [{"Omics": LAB[o], "Entities uploaded": 100, "In the seed list": 40, "Seeds on STRING": 35, "Not on STRING": 5, "Share of the combined weight": 0.25, "Used": True} for o in OM],
    "omics_used": OM, "seeds": {"n_seeds": 14, "n_multi_omics_seeds": 3, "weighting": "equal weight per omics", "n_seeds_not_in_string": 0},
    "metabolites": metab,
    "propagation": {"table": {"columns": cols, "rows": hit_rows}, "n_hits": 22, "n_hits_seed": 14, "n_recruited": 8, "n_recruited_all_family": 8, "n_proteins": 16201,
                    "n_seeds_mapped": 14, "n_seeds_not_in_string": 0,
                    "null_diagnostics": {"n_exceed_top": 0, "null_skew_at_max": 2.1, "n_gamma_tail": 1200, "n_permutations": 500, "skew_caveat": "The null is right-skewed."}, "edges": edges},
    "modules": {"table": {"columns": list(mod_rows[0].keys()), "rows": mod_rows}, "membership": {"M1": [f"9606.P{i}" for i in range(0, 9)], "M2": [f"9606.P{i}" for i in range(10, 16)]},
                "universe_size": 16201, "caveat": None, "selection_caveat": "Labels are descriptive.", "cap_info": {"capped": False, "n_used_for_modules": 22, "n_hits": 22},
                "omics": OM,
                "enrichment": {m: {"namespace": "biological_process", "table": enr_tbl, "candidates": {"biological_process": {"table": enr_tbl}}} for m in ("M1", "M2")}},
    "caveats": ["Equal weight per omics: each omics' seeds are scaled to sum to 1."], "pathway_network": json.loads(json.dumps(pw)),
}
check(len(pw["nodes"]) >= 25 and len(pw["edges"]) > 5 and len(pw["groups"]) >= 2, f"(fixture) the engine's pathway network has {len(pw['nodes'])} nodes, {len(pw['edges'])} links, {len(pw['groups'])} groups")
_hits = {o: sum(1 for n in pw["nodes"] if n["hits"].get(o, 0) > 0) for o in pw["omics"]}
check(pw["n_merged"] >= 1 and any(n["merged"] for n in pw["nodes"]), f"(fixture) near-identical pathways were merged ({pw['n_merged']})")
check(0 < min(_hits.values()) < len(pw["nodes"]) and any(v < len(pw["nodes"]) for v in _hits.values()), f"(fixture) some pathways are hit by an omics and some are not ({_hits} of {len(pw['nodes'])})")

INJECT = """(payload) => {
    STATE.view = 'annotation_analytics';
    STATE.multiomics = freshMultiomicsState();
    STATE.multiomics.method = 'network';
    STATE.multiomics.selected = ['genomics', 'transcriptomics', 'proteomics', 'metabolomics'];
    STATE.multiomics.resolved = {};
    STATE.multiomics.net.run = {results: payload, busy: false};
    STATE.aaPage = 'multiomics_net_results';
    render();
}"""
HAS_DATA = """(id) => { const d = document.getElementById(id); return !!(d && d.data && d.data.length); }"""

with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_context(viewport={"width": 1400, "height": 1000}, accept_downloads=True).new_page()
    errors = []
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    page.on("console", lambda msg: errors.append("console: " + msg.text) if msg.type == "error" else None)
    page.goto(BASE, wait_until="networkidle")
    page.evaluate(INJECT, payload)
    page.wait_for_timeout(1500)
    rail = [t.strip() for t in page.locator("button.result-nav-btn").all_inner_texts()]
    check(rail == ["Overview", "Hit table", "Figures", "Modules", "Module enrichment", "Metabolite enzymes", "Pathway network", "Parameters used"], f"all panes (got {rail})")

    print("=== hit table ===")
    page.locator("button.result-nav-btn", has_text="Hit table").click()
    page.wait_for_timeout(500)
    cb = page.locator("input.multiomics-net-recruited-only")
    check(cb.is_checked() and page.locator(".result-pane:visible tbody tr").count() == 8, "it opens on the RECRUITED view: the 8 non-seed significant proteins")
    cb.click()
    page.wait_for_timeout(400)
    check(page.locator(".result-pane:visible tbody tr").count() == 40, "unticking shows all 40 proteins")
    check(page.locator(".result-pane:visible tbody tr.row-significant").count() == 22, "the 22 significant proteins are highlighted")
    cell = page.locator(".result-pane:visible tbody tr").nth(3).inner_text()
    check("Transcriptomics" in cell and "Proteomics" in cell, "a seed hit by two omics names both")

    print("=== figures ===")
    page.locator("button.result-nav-btn", has_text="Figures").click()
    page.wait_for_timeout(500)
    page.locator("button", has_text="Show network figure").click()
    page.locator("button", has_text="Show null distribution").click()
    page.wait_for_timeout(2000)
    check(page.locator(".result-pane:visible .js-plotly-plot").count() == 2, "the single-omics tiles' network and null-distribution figures draw from this payload")

    print("=== modules ===")
    page.locator("button.result-nav-btn", has_text="Modules").first.click()
    page.wait_for_timeout(900)
    mt = [t.strip() for t in page.locator(".result-pane:visible thead th").all_inner_texts()]
    check(all(f"Seeds from {LAB[o]}" in mt for o in OM), "the module table has a seed column for every omics")
    check(page.evaluate(HAS_DATA, "multiomics-net-bars"), "the composition chart is drawn")
    bars = page.evaluate("""() => { const d = document.getElementById('multiomics-net-bars');
        return {names: d.data.map(t => t.name), x: d.data.map(t => t.x), y: d.data[0].y, stack: d._fullLayout.barmode, tool: !!d.querySelector('.modebar')}; }""")
    check(bars["names"] == [LAB[o] for o in OM] + ["Recruited (no omics supplied them)"] and bars["stack"] == "stack", f"one stacked series per omics plus the recruited share ({bars['names']})")
    check(bars["y"] == ["M1", "M2"] and bars["x"][0] == [2, 0] and bars["x"][1] == [4, 1] and bars["x"][4] == [3, 2], "the bar lengths equal the module table's own counts (M1: 2 genomics, 4 transcript; recruited 3 / 2)")
    check(bars["tool"], "the composition chart has the figure tools row")
    page.locator("button", has_text="Show module-coloured network").click()
    page.wait_for_timeout(1500)
    check(page.locator(".result-pane:visible .js-plotly-plot").count() >= 2, "the module-coloured network draws")

    print("=== module enrichment ===")
    page.locator("button.result-nav-btn", has_text="Module enrichment").click()
    page.wait_for_timeout(700)
    check(page.locator("select.multiomics-net-module-select option").count() == 2 and page.locator(".result-pane:visible tbody tr").count() == 6, "a module selector and the module's 6-term table")
    page.locator("button", has_text="Show dot plot").click()
    page.wait_for_timeout(1500)
    check(page.locator(".result-pane:visible .js-plotly-plot").count() == 1, "the dot plot draws")
    page.locator("select.multiomics-net-module-select").select_option("M2")
    page.wait_for_timeout(500)
    check(page.locator(".result-pane:visible tbody tr").count() == 6, "switching module redraws the table")

    print("=== metabolites ===")
    page.locator("button.result-nav-btn", has_text="Metabolite enzymes").click()
    page.wait_for_timeout(500)
    txt = page.inner_text(".result-pane:visible")
    check("left out: in 2103 reactions" in txt and page.locator(".result-pane:visible tbody tr").count() == 2, "the hub is listed with its reason")

    print("=== pathway network ===")
    page.locator("button.result-nav-btn", has_text="Pathway network").click()
    page.wait_for_timeout(1800)
    m = page.evaluate("""(nodes) => { const d = document.getElementById('multiomics-net-pathway'); const xa = d._fullLayout.xaxis, ya = d._fullLayout.yaxis;
        const nodeT = d.data.find(t => t.mode === 'markers+text');
        const pos = nodeT.x.map((x, i) => [x, nodeT.y[i]]);
        const labels = [...d.querySelectorAll('.textpoint text')].map(t => { const r = t.getBoundingClientRect(); return [r.left, r.top, r.right, r.bottom]; });
        let overlaps = 0;
        for (let i = 0; i < labels.length; i++) for (let j = i + 1; j < labels.length; j++) { const a = labels[i], b = labels[j];
            if (a[0] < b[2] - 1 && b[0] < a[2] - 1 && a[1] < b[3] - 1 && b[1] < a[3] - 1) overlaps++; }
        const chipSeries = d.data.filter(t => t.mode === 'markers' && t.name).map(t => [t.name, t.x.length]);
        const chipPos = Object.fromEntries(d.data.filter(t => t.mode === 'markers' && t.name).map(t => [t.name, [t.x, t.y]]));
        const gT = d.data.find(t => t.name === 'Genomics' && t.mode === 'markers');
        const withG = nodes.map((n, i) => [n, i]).filter(([n]) => (n.hits.genomics || 0) > 0);
        const belowPx = withG.map(([n, i], k) => ya.l2p(gT.y[k]) - ya.l2p(nodeT.y[i]));
        const outlines = nodeT.marker.line.width.filter(w => w > 2).length;
        const colours = new Set(nodeT.marker.color).size;
        return {n: nodeT.x.length, pos, chipPos, belowPx, overlaps, nLabels: labels.length, chipSeries, outlines, colours, tool: !!d.querySelector('.modebar'), w: d.clientWidth, h: d.clientHeight,
                scroll: getComputedStyle(d.parentElement).overflow, edges: d.data.filter(t => t.mode === 'lines').length}; }""", payload["pathway_network"]["nodes"])
    nodes = payload["pathway_network"]["nodes"]
    check(m["n"] == len(nodes) and m["nLabels"] == len(nodes), f"all {len(nodes)} pathways are drawn with a label ({m['n']} / {m['nLabels']})")
    check(m["w"] == payload["pathway_network"]["canvas"]["width"], f"the chart is exactly the width the engine sized the canvas to ({m['w']}px)")
    check(all(abs(a[0] - n["x"]) < 1e-9 and abs(a[1] - n["y"]) < 1e-9 for a, n in zip(m["pos"], nodes)), "the nodes sit at the engine's own layout positions")
    hov = page.evaluate("document.getElementById('multiomics-net-pathway').data.find(t => t.mode === 'markers+text').hovertext")
    merged_nodes = [i for i, n in enumerate(nodes) if n["merged"]]
    check(all("also stands for" in hov[i] for i in merged_nodes) and not any("also stands for" in hov[i] for i, n in enumerate(nodes) if not n["merged"]),
          f"SENSITIVITY: the hover of a pathway that absorbed others says so ({len(merged_nodes)} of {len(nodes)} pathways)")
    check(f"{payload['pathway_network']['n_merged']} further pathway(s) were merged" in page.inner_text(".result-pane:visible"), "the pane says how many pathways were merged")
    check(m["overlaps"] == 0, f"NO TWO LABELS OVERLAP on a crowded network of long pathway names (measured: {m['overlaps']} overlapping pairs)")
    shown = page.evaluate("[...document.querySelectorAll('#multiomics-net-pathway .textpoint text')].map(t => t.textContent).join(' ')")
    check("<br>" not in shown and "&lt;" not in shown and "&amp;" not in shown, "SENSITIVITY: the labels show no visible '<br>' or HTML entity (a wrapped name is real lines)")
    lines = page.evaluate("Math.max(...[...document.querySelectorAll('#multiomics-net-pathway .textpoint text')].map(t => t.querySelectorAll('tspan').length))")
    check(lines >= 2, f"long names are wrapped onto several lines (up to {lines} lines)")
    have = {LAB[o]: sum(1 for n in nodes if n["hits"].get(o, 0) > 0) for o in payload["pathway_network"]["omics"]}
    check(dict(m["chipSeries"]) == have, f"chips appear only where an omics hit the pathway ({have})")
    check(m["belowPx"] and min(m["belowPx"]) > 8, f"SENSITIVITY: AS DRAWN, the chips are under their node on the page (at least {min(m['belowPx']):.0f}px lower)")
    gen_nodes = [n for n in nodes if n["hits"].get("genomics", 0) > 0]
    gx, gy = m["chipPos"]["Genomics"]
    dys = [cy - n["y"] for cy, n in zip(gy, gen_nodes)]
    dxs = [cx - n["x"] for cx, n in zip(gx, gen_nodes)]
    check(all(10 < d < 45 for d in dys) and all(abs(d) < 30 for d in dxs), f"SENSITIVITY: the chips sit just BELOW their node (offsets {min(dys):.0f}..{max(dys):.0f}px down, within {max(abs(d) for d in dxs):.0f}px sideways)")
    check(m["outlines"] == sum(1 for n in nodes if n["significant"]), f"significant pathways are outlined ({m['outlines']})")
    check(m["colours"] == len(payload["pathway_network"]["groups"]) or m["colours"] >= 2, f"nodes are coloured by group ({m['colours']} colours, {len(payload['pathway_network']['groups'])} groups)")
    check(m["edges"] == len(payload["pathway_network"]["edges"]), f"one line per link ({m['edges']})")
    check(m["tool"] and m["scroll"] == "auto" and m["w"] >= 1100, f"figure tools row; fixed size {m['w']}x{m['h']} inside a scroll frame")
    page.screenshot(path=os.path.join(os.environ.get("NET_SHOT_DIR", _ROOT), "net_chart.png")) if os.environ.get("NET_SHOT_DIR") else None

    print("=== the figure tools: the same options as every other figure ===")
    W0 = payload["pathway_network"]["canvas"]["width"]
    H0 = payload["pathway_network"]["canvas"]["height"]
    RANGE = "() => { const d = document.getElementById('multiomics-net-pathway'); return {x: d._fullLayout.xaxis.range.slice(), y: d._fullLayout.yaxis.range.slice(), drag: d._fullLayout.dragmode, w: d.clientWidth}; }"
    tool_names = ["Zoom tool", "Pan tool", "Zoom in", "Zoom out", "Autoscale", "Reset axes"]
    btns = page.locator(".result-pane:visible .multiomics-net-path-tools button")
    check([t.strip() for t in btns.all_inner_texts()] == tool_names, f"a visible tools row has {tool_names}")
    vp = page.evaluate("window.innerWidth")
    in_view = page.evaluate("""() => [...document.querySelectorAll('.result-pane .multiomics-net-path-tools button, .multiomics-net-path-download, select.multiomics-net-path-size, select.multiomics-net-path-labels')]
        .filter(e => e.offsetParent !== null).map(e => { const r = e.getBoundingClientRect(); return [r.left, r.right]; })""")
    check(len(in_view) == 9 and all(l >= 0 and r <= vp for l, r in in_view), f"SENSITIVITY: every tool, the Size and name menus and the Download button are on screen without scrolling sideways (the chart's own toolbar is {W0}px away)")
    r0 = page.evaluate(RANGE)
    page.locator(".result-pane:visible button.fig-tool-in").first.click()
    page.wait_for_timeout(500)
    r1 = page.evaluate(RANGE)
    check(abs((r1["x"][1] - r1["x"][0]) - (r0["x"][1] - r0["x"][0]) / 2) < 1, f"Zoom in halves the visible width ({r0['x'][1] - r0['x'][0]:.0f} -> {r1['x'][1] - r1['x'][0]:.0f} units)")
    page.locator(".result-pane:visible button.fig-tool-out").first.click()
    page.wait_for_timeout(500)
    r2 = page.evaluate(RANGE)
    check(abs((r2["x"][1] - r2["x"][0]) - (r0["x"][1] - r0["x"][0])) < 1, "Zoom out puts the width back")
    page.locator(".result-pane:visible button.fig-tool-zoom").first.click()
    page.wait_for_timeout(300)
    check(page.evaluate(RANGE)["drag"] == "zoom", "Zoom tool switches dragging to box-zoom")
    page.locator(".result-pane:visible button.fig-tool-pan").first.click()
    page.wait_for_timeout(300)
    check(page.evaluate(RANGE)["drag"] == "pan", "Pan tool switches dragging to move")
    page.locator(".result-pane:visible button.fig-tool-in").first.click()
    page.locator(".result-pane:visible button.fig-tool-in").first.click()
    page.wait_for_timeout(500)
    page.locator(".result-pane:visible button.fig-tool-reset").first.click()
    page.wait_for_timeout(500)
    r3 = page.evaluate(RANGE)
    check(r3["x"] == [0, W0] and r3["y"] == [H0, 0], f"Reset axes returns to the opening view ({r3['x']}, {r3['y']})")
    page.locator(".result-pane:visible button.fig-tool-in").first.click()
    page.wait_for_timeout(400)
    page.locator(".result-pane:visible button.fig-tool-auto").first.click()
    page.wait_for_timeout(500)
    r4 = page.evaluate(RANGE)
    check(r4["x"] != r3["x"] or r4["y"] != r3["y"], "Autoscale fits the axes to the data (a different view from the zoomed one)")

    SIZES = """() => { const d = document.getElementById('multiomics-net-pathway'); const xa = d._fullLayout.xaxis, ya = d._fullLayout.yaxis;
        const nodeT = d.data.find(t => t.mode === 'markers+text'); const gT = d.data.find(t => t.name === 'Genomics' && t.mode === 'markers');
        const nodes = %s; const withG = nodes.map((n, i) => [n, i]).filter(([n]) => (n.hits.genomics || 0) > 0);
        const below = withG.map(([n, i], k) => ya.l2p(gT.y[k]) - ya.l2p(nodeT.y[i]));
        return {w: d.clientWidth, h: d.clientHeight, below: [Math.min(...below), Math.max(...below)], frame: d.parentElement.clientWidth,
                labels: [...d.querySelectorAll('.textpoint text')].filter(t => t.textContent.trim()).length}; }""" % json.dumps(nodes)
    for val, want in (("0.5", 0.5), ("1.5", 1.5), ("2", 2.0), ("1", 1.0)):
        page.locator("select.multiomics-net-path-size").select_option(val)
        page.wait_for_timeout(1500)
        sz = page.evaluate(SIZES)
        check(abs(sz["w"] - round(W0 * want)) <= 2, f"Size {int(want * 100)}%: the chart is {sz['w']}px wide (canvas {W0}px x {want})")
        check(10 < sz["below"][0] and sz["below"][1] < 45, f"...and the chips still sit just under their node at that size ({sz['below'][0]:.0f}..{sz['below'][1]:.0f}px)")
    page.locator("select.multiomics-net-path-size").select_option("fit")
    page.wait_for_timeout(1500)
    sz = page.evaluate(SIZES)
    check(sz["w"] <= sz["frame"] and sz["w"] > 300, f"Size: fit to page makes the chart fit its frame ({sz['w']}px in {sz['frame']}px)")
    page.locator("select.multiomics-net-path-size").select_option("1")
    page.wait_for_timeout(1200)
    page.locator("select.multiomics-net-path-labels").select_option("none")
    page.wait_for_timeout(1200)
    check(page.evaluate(SIZES)["labels"] == 0, "Pathway names: none - no name is drawn")
    page.locator("select.multiomics-net-path-labels").select_option("top")
    page.wait_for_timeout(1200)
    check(page.evaluate(SIZES)["labels"] == 10, "Pathway names: the 10 best only - exactly ten names")
    page.locator("select.multiomics-net-path-labels").select_option("all")
    page.wait_for_timeout(1200)
    check(page.evaluate(SIZES)["labels"] == len(nodes), f"Pathway names: all - all {len(nodes)}")
    check(page.evaluate(HAS_DATA, "multiomics-net-pathway") and page.evaluate("document.getElementById('multiomics-net-pathway')._figHome.x[1]") == W0, "the redrawn chart keeps its home view")
    page.locator(".multiomics-net-path-download").click()
    page.wait_for_timeout(500)
    menu = page.evaluate("""() => { const m = document.querySelector('.ora-chart-export-menu'); const b = document.querySelector('.multiomics-net-path-download');
        if (!m) return null; const mr = m.getBoundingClientRect(), br = b.getBoundingClientRect(); return {top: mr.top, btnBottom: br.bottom, left: mr.left, btnLeft: br.left, items: m.innerText}; }""")
    check(menu is not None and abs(menu["top"] - menu["btnBottom"]) < 12 and abs(menu["left"] - menu["btnLeft"]) < 12, f"Download figure opens the app's export menu under the button ({menu['items'].replace(chr(10), ' | ')[:80] if menu else None})")
    check(menu is not None and "PNG" in menu["items"].upper(), "...offering the usual image formats")

    check(len(errors) == 0, f"no page error or console error (got {errors})")
    browser.close()

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
