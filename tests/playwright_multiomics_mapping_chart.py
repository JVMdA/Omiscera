"""Fast browser test of the Cross-omics Pathway Mapping CHART (static/app.js) on SYNTHETIC data.

The server is needed only to serve app.js (no network, no analysis): a hand-built mapping result is put into
STATE and the real render path draws it, so every number below is known in advance and the checks can fail for
the reason they name. playwright_multiomics_mapping.py drives the same chart end to end on real data; this file
exists because the chart's hardest promise - evidence chips that never cover each other - deserves checks that
run in seconds and can be mutation-tested.

What is asserted:
  - the shape of every chip: a triangle for an expression-omics hit (up / down), a diamond when the box's
    members disagree, a filled circle for a list hit, a HOLLOW circle (its own trace) for measured-but-not-
    significant, and no chip at all for an absent state;
  - chips sit in FIXED slots (Genomics, Transcriptomics, Proteomics, Metabolomics from left to right) so an
    omics that is absent on a box leaves its slot empty instead of shifting the others;
  - no two chips on the same box are closer than a marker width, at two window widths, zoomed in 5x, and zoomed
    back out (the positions are in map units, the spacing must stay in screen pixels);
  - chips sit BELOW their box; a box with no evidence has no label; a hit box has a darker border;
  - reaction lines: KEGG-style (substrate to product between boxes) and Reactome-style (the diagram's own segments
    plus input / output / catalyst links to the reaction point, catalysts dotted);
  - an empty result says so instead of drawing nothing.
"""
import os

from playwright.sync_api import sync_playwright

BASE = "http://127.0.0.1:8009"
PASS = FAIL = 0
MARKER = 9


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def member(label, hit=True, log2fc=None, padj=None, via="UniProt"):
    return {"label": label, "hit": hit, "log2fc": log2fc, "padj": padj, "via": via}


def absent():
    return {"state": "absent", "direction": None, "members": []}


def status(**kw):
    out = {o: absent() for o in ("genomics", "transcriptomics", "proteomics", "metabolomics")}
    out.update(kw)
    return out


def node(nid, kind, label, x, y, st, w=40, h=17):
    return {"id": nid, "kind": kind, "label": label, "x": x, "y": y, "w": w, "h": h, "status": st}


OMICS = ["genomics", "transcriptomics", "proteomics", "metabolomics"]
NODES = [
    # three chips: list hit, expression UP, expression DOWN
    node("k1", "gene", "AAA", 100, 100, status(
        genomics={"state": "hit", "direction": None, "members": [member("AAA")]},
        transcriptomics={"state": "hit", "direction": "up", "members": [member("AAA", True, 2.5, 0.001)]},
        proteomics={"state": "hit", "direction": "down", "members": [member("AAA", True, -2.0, 0.003)]})),
    # Transcriptomics ABSENT: Genomics and Proteomics must keep their own slots (not shift together)
    node("k2", "gene", "BBB", 200, 100, status(
        genomics={"state": "hit", "direction": None, "members": [member("BBB")]},
        proteomics={"state": "hit", "direction": "mixed", "members": [member("B1", True, 2.0, 0.01), member("B2", True, -2.0, 0.01)]})),
    # measured, not significant
    node("k3", "gene", "CCC", 300, 100, status(
        genomics={"state": "hit", "direction": None, "members": [member("CCC")]},
        transcriptomics={"state": "measured", "direction": None, "members": [member("CCC", False, 0.2, 0.7)]},
        proteomics={"state": "measured", "direction": None, "members": [member("CCC", False, -0.1, 0.9)]})),
    # a compound hit by Metabolomics only
    node("k4", "compound", "Citrate", 100, 200, status(
        metabolomics={"state": "hit", "direction": None, "members": [member("Citrate", via="KEGG compound ID")]}), w=8, h=8),
    # no evidence at all
    node("k5", "gene", "DDD", 200, 200, status()),
    # a link to another pathway
    node("k6", "link", "Glycolysis", 300, 200, status(), w=60, h=20),
    # four chips on one small compound box: every slot used
    node("k7", "compound", "Pyruvate", 200, 300, status(
        genomics={"state": "hit", "direction": None, "members": [member("x")]},
        transcriptomics={"state": "hit", "direction": "up", "members": [member("x", True, 1.5, 0.01)]},
        proteomics={"state": "hit", "direction": "down", "members": [member("x", True, -1.5, 0.01)]},
        metabolomics={"state": "hit", "direction": None, "members": [member("Pyruvate")]}), w=8, h=8),
]
KEGG_RESULT = {
    "status": "ok", "database": "KEGG", "pathway_id": "hsa99999", "pathway_name": "Synthetic KEGG", "omics": OMICS,
    "nodes": NODES, "note": None, "restricted_to": None, "caveat": "Pathway mapping is descriptive: it computes no p-value.",
    "reactions": [{"id": "rx0", "name": "R1", "reversible": False, "enzymes": ["k1"], "substrates": ["k4"], "products": ["k5"],
                   "pos": None, "segments": []}],
    "summary": {"omics": OMICS, "n_nodes": 7, "n_nodes_with_hit": 5, "n_links": 1,
                "per_omics": {o: {"label": o, "nodes_hit": 1, "nodes_measured_only": 0, "members_hit": 1} for o in OMICS}},
    "members": {"columns": ["Pathway node", "Node type", "Member", "Omics", "Status", "Direction", "log2FC", "padj", "Matched by"], "rows": []},
    "links": {"columns": ["Reaction"], "rows": []},
}
REACTOME_RESULT = dict(KEGG_RESULT, database="Reactome", pathway_id="R-HSA-1", pathway_name="Synthetic Reactome", reactions=[
    {"id": "rx1", "name": "R2", "reversible": False, "enzymes": ["k1"], "substrates": ["k4"], "products": ["k5"],
     "pos": [150, 150], "segments": [[140, 150, 150, 150], [150, 150, 160, 150]]}])



def hit(label, omics_dir=None, via="KEGG gene ID"):
    return {"state": "hit", "direction": omics_dir, "members": [member(label, True, 2.0 if omics_dir == "up" else (-2.0 if omics_dir == "down" else None), 0.01 if omics_dir else None, via)]}


def measured(label):
    return {"state": "measured", "direction": None, "members": [member(label, False, 0.2, 0.7)]}


def tail(result_status, **kw):
    out = {"status": "ok", "database": "KEGG", "pathway_id": "hsa00001", "pathway_name": "Synthetic", "omics": OMICS, "note": None,
           "restricted_to": None, "caveat": "Pathway mapping is descriptive: it computes no p-value.", "reactions": [], "relations": [],
           "decor_lines": [], "line_graphics": False, "interactions": None, "subpathways": None,
           "summary": {"omics": OMICS, "n_nodes": 0, "n_nodes_with_hit": 0, "n_links": 0, "n_interactions": 0,
                       "per_omics": {o: {"label": o, "nodes_hit": 0, "nodes_measured_only": 0, "members_hit": 0} for o in OMICS}},
           "members": {"columns": ["Pathway node", "Box ID", "Node type", "Member", "Omics", "Status", "Direction", "log2FC", "padj", "Matched by"], "rows": []},
           "links": {"columns": ["Reaction"], "rows": []}}
    out.update(kw)
    return out


# KEGG overview map: genes are LINES (a polyline, flat x1, y1, x2, y2, ...), coloured by their pathway CLASS
CARB, LIPID = "Carbohydrate metabolism", "Lipid metabolism"
LINE_NODES = [
    node("k1", "gene", "G1only", 200, 0, status(genomics=hit("G1only")), w=0, h=0),                                   # carb, hit
    node("k2", "gene", "G12", 200, 100, status(genomics=hit("G12"), transcriptomics=hit("G12", "up")), w=0, h=0),      # lipid, hit by two
    node("k3", "gene", "Gmeas", 200, 200, status(transcriptomics=measured("Gmeas")), w=0, h=0),                        # carb, measured only
    node("k4", "gene", "Gnone", 200, 300, status(), w=0, h=0),                                                         # lipid, untouched
    node("k5", "compound", "Cmp", 0, 0, status(metabolomics=hit("Cmp", None, "KEGG compound ID")), w=8, h=8),
    node("k6", "compound", "Cquiet", 400, 300, status(), w=8, h=8),
    node("k7", "link", "Glycolysis", 300, 350, status(), w=60, h=20),
]
for nid, y, cls in (("k1", 0, CARB), ("k2", 100, LIPID), ("k3", 200, CARB), ("k4", 300, LIPID)):
    nn = next(n for n in LINE_NODES if n["id"] == nid)
    nn["path"] = [0, y, 200, y, 400, y]
    nn["cls"] = cls
next(n for n in LINE_NODES if n["id"] == "k7")["cls"] = CARB
LINE_RESULT = tail("ok", nodes=LINE_NODES, line_graphics=True, decor_lines=[[0, 350, 400, 350]],
                   classes=[{"name": CARB, "color": "#B3B3E6", "n_lines": 2, "n_hit": 1}, {"name": LIPID, "color": "#80CCCC", "n_lines": 2, "n_hit": 1}],
                   reactions=[{"id": "rx0", "name": "R1", "reversible": False, "enzymes": ["k1"], "substrates": ["k5"], "products": ["k6"], "pos": None, "segments": []}])
WIDE_RESULT = tail("ok", nodes=[node("w1", "gene", "W1", 0, 0, status(), w=40, h=17), node("w2", "gene", "W2", 2000, 1200, status(genomics=hit("W2")), w=40, h=17)])

# KEGG signalling map: relations (arrow / bar / plain / dashed ...), a group, a gene hit, one not hit
SIG_NODES = [
    node("a", "gene", "AAA", 100, 100, status(genomics=hit("AAA"), transcriptomics=hit("AAA", "up"))),
    node("b", "gene", "BBB", 300, 100, status(genomics=hit("BBB"))),
    node("c", "gene", "CCC", 500, 100, status()),
    node("e", "gene", "EEE", 600, 300, status(proteomics=hit("EEE", "down"))),
    node("g", "group", "CCC+EEE", 550, 200, status(proteomics=hit("EEE", "down")), w=130, h=30),
    node("d", "compound", "Cit", 100, 300, status(metabolomics=hit("Cit", None, "KEGG compound ID")), w=8, h=8),
]
SIG_NODES[4]["components"] = ["c", "e"]
SIG_RESULT = tail("ok", nodes=SIG_NODES, relations=[
    {"a": "a", "b": "b", "rtype": "pprel", "types": ["activation"]},                                  # both hit: dark arrow
    {"a": "b", "b": "c", "rtype": "pprel", "types": ["inhibition"]},                                  # CCC not hit: light bar
    {"a": "a", "b": "g", "rtype": "pprel", "types": ["binding/association"]},                         # group hit: dark plain
    {"a": "b", "b": "e", "rtype": "pprel", "types": ["activation", "phosphorylation"]},               # dark arrow + '+p'
    {"a": "c", "b": "e", "rtype": "pprel", "types": ["indirect effect"]},                              # light dashed arrow
    {"a": "e", "b": "a", "rtype": "pprel", "types": ["dissociation"]},                                 # dark dash-dot
    {"a": "a", "b": "nowhere", "rtype": "pprel", "types": ["activation"]},                            # skipped
], interactions={"columns": ["From", "To", "Relation", "From hit by", "To hit by"],
                 "rows": [{"From": "AAA", "To": "BBB", "Relation": "activation", "From hit by": "Genomics; Transcriptomics (up)", "To hit by": "Genomics"}]})

# a crowded map: over 300 boxes, so the default is to label only the hits
BIG_NODES = [node(f"n{i}", "gene", f"Gene{i}", 50 * (i % 20), 40 * (i // 20), status(genomics=hit(f"Gene{i}")) if i % 50 == 0 else
                  (status(transcriptomics=measured(f"Gene{i}")) if i % 50 == 1 else status())) for i in range(320)]
BIG_RESULT = tail("ok", nodes=BIG_NODES)

CHIP_MIN_JS = """() => {
    const d = document.getElementById('multiomics-pathway-map'); const xa = d._fullLayout.xaxis, ya = d._fullLayout.yaxis;
    const byNode = {};
    d.data.forEach((t) => { if (!t.customdata || t.mode !== 'markers' || !t.x) return;
        t.customdata.forEach((id, i) => { if (typeof id === 'string' && t.x[i] !== null) (byNode[id] = byNode[id] || []).push([xa.l2p(t.x[i]), ya.l2p(t.y[i])]); }); });
    let minDist = 1e9, three = 0;
    Object.values(byNode).forEach(pts => { if (pts.length >= 3) three++;
        for (let a = 0; a < pts.length; a++) for (let b = a + 1; b < pts.length; b++) minDist = Math.min(minDist, Math.hypot(pts[a][0] - pts[b][0], pts[a][1] - pts[b][1])); });
    return {minDist, three};
}"""
PIXELS_JS = """(id) => {
    const d = document.getElementById('multiomics-pathway-map'); const fl = d._fullLayout, xa = fl.xaxis, ya = fl.yaxis;
    const rect = d.getBoundingClientRect(); const p = d._mapPos[id];
    return {x: rect.left + xa._offset + xa.l2p(p.x), y: rect.top + ya._offset + ya.l2p(p.y), px: p.x, py: p.y};
}"""
SNAP_JS = """() => {
    const d = document.getElementById('multiomics-pathway-map'); const fl = d._fullLayout, xa = fl.xaxis, ya = fl.yaxis;
    const P = (x, y) => [xa._offset + xa.l2p(x), ya._offset + ya.l2p(y)];
    const traces = d.data.map((t, i) => ({i, mode: t.mode, name: t.name || '', dash: t.line && t.line.dash, color: t.line && t.line.color, width: t.line && t.line.width,
        x: (t.x || []).slice(), y: (t.y || []).slice(), text: t.text, sym: t.marker && t.marker.symbol, msize: t.marker && t.marker.size, custom: t.customdata,
        mcolor: t.marker && t.marker.color, hover: t.hovertext, tcolor: t.textfont && t.textfont.color}));
    return {traces, ppu: Math.abs(xa.l2p(1) - xa.l2p(0)), scrollZoom: d._context && d._context.scrollZoom, dragmode: d._fullLayout.dragmode,
            pos: JSON.parse(JSON.stringify(d._mapPos))};
}"""

METRICS_JS = """() => {
    const d = document.getElementById('multiomics-pathway-map');
    if (!d || !d._fullLayout || !d.data) return null;
    const xa = d._fullLayout.xaxis, ya = d._fullLayout.yaxis;
    const byNode = {}, nodeTraces = [], chipTraces = [];
    d.data.forEach((t, ti) => {
        if (t.mode === 'markers+text') nodeTraces.push(ti);
        if (!t.customdata || t.mode !== 'markers' || !t.x) return;
        chipTraces.push(ti);
        const sym = t.marker.symbol;
        t.customdata.forEach((id, i) => {
            if (typeof id !== 'string' || t.x[i] === null) return;
            (byNode[id] = byNode[id] || []).push({name: t.name || '(measured)', x: xa.l2p(t.x[i]), y: ya.l2p(t.y[i]),
                symbol: Array.isArray(sym) ? sym[i] : sym, trace: ti});
        });
    });
    let minDist = 1e9;
    Object.values(byNode).forEach(pts => { for (let a = 0; a < pts.length; a++) for (let b = a + 1; b < pts.length; b++)
        minDist = Math.min(minDist, Math.hypot(pts[a].x - pts[b].x, pts[a].y - pts[b].y)); });
    const nodePx = {};
    nodeTraces.forEach(ti => { const t = d.data[ti]; t.hovertext.forEach((h, i) => {
        const m = /<b>(.*?)<\\/b>/.exec(h); if (m) nodePx[m[1]] = {x: xa.l2p(t.x[i]), y: ya.l2p(t.y[i]), text: t.text[i], i, ti,
            border: t.marker.line.color[i], width: t.marker.line.width[i]}; }); });
    return {minDist, byNode, nodePx, ppu: Math.abs(xa.l2p(1) - xa.l2p(0)), width: d.clientWidth,
            names: d.data.map(t => t.name || ''), lines: d.data.filter(t => t.mode === 'lines').map(t => ({n: t.x.length, dash: t.line.dash || 'solid'}))};
}"""


def show(page, result):
    page.evaluate("""(res) => {
        STATE.view = 'annotation_analytics'; STATE.aaPage = 'multiomics_map'; STATE._resultNavActive = {};
        const mo = freshMultiomicsState(); mo.method = 'mapping'; mo.selected = ['genomics','transcriptomics','proteomics','metabolomics'];
        mo.resolved = {omics: {}}; mo.taxonId = '9606'; mo.taxonName = 'Homo sapiens'; mo.map.result = res; mo.map.database = res.database;
        STATE.multiomics = mo; render(); }""", result)
    page.wait_for_timeout(1500)


with sync_playwright() as p:
    browser = p.chromium.launch()
    for width in (1400, 800):
        context = browser.new_context(viewport={"width": width, "height": 1000})
        page = context.new_page()
        errors = []
        page.on("pageerror", lambda exc: errors.append(str(exc)))
        page.goto(BASE, wait_until="networkidle")
        page.wait_for_timeout(300)
        print(f"=== KEGG-style synthetic map, window {width}px ===")
        show(page, KEGG_RESULT)
        m = page.evaluate(METRICS_JS)
        check(m is not None, f"[{width}] the chart is drawn from a mapping result")
        by = m["byNode"]
        check(set(by) == {"k1", "k2", "k3", "k4", "k7"}, f"[{width}] chips only on boxes with evidence ({sorted(by)}); none on the empty box or the link")
        sym = {k: {c["name"]: c["symbol"] for c in v} for k, v in by.items()}
        check(sym["k1"] == {"Genomics": "circle", "Transcriptomics": "triangle-up", "Proteomics": "triangle-down"},
              f"[{width}] list hit = circle, expression UP = triangle-up, DOWN = triangle-down ({sym['k1']})")
        check(sym["k2"].get("Proteomics") == "diamond", f"[{width}] a box whose members disagree is a diamond ({sym['k2']})")
        check(sym["k3"].get("(measured)") == "circle-open" and "Transcriptomics" not in sym["k3"],
              f"[{width}] SENSITIVITY: measured-not-significant is a HOLLOW circle in its own trace, not a filled hit ({sym['k3']})")
        check(sorted(sym["k7"]) == ["Genomics", "Metabolomics", "Proteomics", "Transcriptomics"], f"[{width}] a box can carry all four chips")
        x_of = lambda nid, name: next(c["x"] for c in by[nid] if c["name"] == name)
        check(x_of("k1", "Genomics") < x_of("k1", "Transcriptomics") < x_of("k1", "Proteomics"),
              f"[{width}] chips run Genomics, Transcriptomics, Proteomics from left to right")
        check(x_of("k7", "Proteomics") < x_of("k7", "Metabolomics") and x_of("k7", "Transcriptomics") < x_of("k7", "Proteomics"),
              f"[{width}] ...and Metabolomics is the fourth slot")
        gap_k1 = x_of("k1", "Proteomics") - x_of("k1", "Genomics")
        gap_k2 = x_of("k2", "Proteomics") - x_of("k2", "Genomics")
        check(abs(gap_k2 - gap_k1) < 1.5 and gap_k2 >= 2 * MARKER,
              f"[{width}] SENSITIVITY: with Transcriptomics ABSENT on a box the other chips keep their slots (gap {gap_k2:.1f}px = the gap with it present, {gap_k1:.1f}px)")
        check(m["minDist"] >= MARKER, f"[{width}] SENSITIVITY: no two chips on one box are closer than a marker ({m['minDist']:.1f}px; stacked would be 0)")
        below = all(c["y"] > m["nodePx"][lbl]["y"] for lbl, nid in (("AAA", "k1"), ("BBB", "k2"), ("CCC", "k3"), ("Citrate", "k4"), ("Pyruvate", "k7"))
                    for c in by[nid])
        check(below, f"[{width}] every chip sits BELOW its box")
        check(m["nodePx"]["DDD"]["text"] == "" and m["nodePx"]["AAA"]["text"] == "AAA",
              f"[{width}] a box with no evidence is not labelled; one with evidence is")
        check(m["nodePx"]["AAA"]["width"] > m["nodePx"]["DDD"]["width"] and m["nodePx"]["AAA"]["border"] != m["nodePx"]["DDD"]["border"],
              f"[{width}] a hit box has a darker, thicker border than an untouched one")
        # zoom in 5x on the cluster, then back out: the pixel spacing must hold at every zoom
        page.evaluate("""() => { const d = document.getElementById('multiomics-pathway-map'); const r = d._fullLayout.xaxis.range;
            const mid = (r[0] + r[1]) / 2, w = (r[1] - r[0]) / 10; Plotly.relayout(d, {'xaxis.range': [mid - w, mid + w]}); }""")
        page.wait_for_timeout(1200)
        mz = page.evaluate(METRICS_JS)
        check(mz["ppu"] > m["ppu"] * 3 and mz["minDist"] >= MARKER,
              f"[{width}] SENSITIVITY: zoomed in ({m['ppu']:.2f} -> {mz['ppu']:.2f} px per map unit) chips stay a marker apart ({mz['minDist']:.1f}px)")
        page.evaluate("""() => { const d = document.getElementById('multiomics-pathway-map');
            Plotly.relayout(d, {'xaxis.autorange': true, 'yaxis.autorange': true}); }""")
        page.wait_for_timeout(1200)
        mo = page.evaluate(METRICS_JS)
        check(mo["minDist"] >= MARKER, f"[{width}] ...and zoomed back out ({mo['minDist']:.1f}px)")
        # zoom OUT far beyond the map: the chips must not collapse onto each other
        page.evaluate("""() => { const d = document.getElementById('multiomics-pathway-map'); const r = d._fullLayout.xaxis.range;
            const mid = (r[0] + r[1]) / 2, w = (r[1] - r[0]) * 2; Plotly.relayout(d, {'xaxis.range': [mid - w, mid + w]}); }""")
        page.wait_for_timeout(1200)
        mf = page.evaluate(METRICS_JS)
        check(mf["ppu"] < m["ppu"] and mf["minDist"] >= MARKER, f"[{width}] SENSITIVITY: zoomed far OUT ({mf['ppu']:.2f} px/unit) chips are still a marker apart ({mf['minDist']:.1f}px)")
        lines = m["lines"]
        check(len(lines) == 2 and lines[0]["n"] == 3 and lines[1]["n"] == 0,
              f"[{width}] KEGG: one substrate-to-product line (3 points incl. the break), no catalyst lines ({lines})")
        check("Glycolysis" not in [n for n in m["nodePx"]] and "Glycolysis" not in page.evaluate("document.getElementById('multiomics-pathway-map').data.filter(t => t.mode === 'markers+text').map(t => t.text.join('|')).join('|')"),
              f"[{width}] a link to another pathway is plain text, never a box")
        check(len(errors) == 0, f"[{width}] no page error (got {errors})")

        if width == 1400:
            print("=== Reactome-style synthetic map ===")
            show(page, REACTOME_RESULT)
            mr = page.evaluate(METRICS_JS)
            lr = mr["lines"]
            check(lr[0]["n"] == 3 * (2 + 2) and lr[1]["n"] == 3, f"Reactome: the diagram's own 2 segments plus input / output links, and a catalyst link ({lr})")
            check(lr[1]["dash"] == "dot", "the catalyst link is dotted")
            check(mr["minDist"] >= MARKER, f"chips clear on the Reactome-style map too ({mr['minDist']:.1f}px)")
            print("=== the result card ===")
            panes = page.locator("button.result-nav-btn").all_inner_texts()
            check(panes == ["Pathway map", "Members", "Reaction links", "Overview"], f"panes: {panes}")
            body = page.inner_text("body")
            check("Synthetic Reactome (R-HSA-1, Reactome)" in body, "headed with the pathway name, ID and database")
            check("filled = hit, hollow = measured but not significant" in body, "the map pane explains the chips")
            page.locator("button.result-nav-btn", has_text="Members").click()
            page.wait_for_timeout(300)
            check("None of your omics' members could be placed" in page.inner_text("body"), "an empty member table says so")
            page.locator("button.result-nav-btn", has_text="Reaction links").click()
            page.wait_for_timeout(300)
            check("No reaction in this pathway has cross-omics evidence" in page.inner_text("body"), "an empty link table says so")
            page.locator("button.result-nav-btn", has_text="Overview").click()
            page.wait_for_timeout(300)
            body = page.inner_text("body")
            check("computes no p-value" in body and "Boxes hit by any omics" in body, "overview: caveat and metric cards")
            print("=== a pathway that could not be fetched / an empty one ===")
            show(page, dict(KEGG_RESULT, status="nodata", note="KEGG has no pathway map (KGML) for hsa99999.", nodes=[]))
            check("KEGG has no pathway map" in page.inner_text("body") and page.locator("#multiomics-pathway-map").count() == 0,
                  "a pathway with no map shows the reason and draws no chart")
            show(page, dict(KEGG_RESULT, nodes=[], reactions=[]))
            check("No map available for this pathway" in page.inner_text("#multiomics-pathway-map"), "an ok result with no boxes says so instead of drawing nothing")
        if width == 1400:

            print("=== KEGG overview map: genes are lines, coloured by pathway class ===")
            show(page, LINE_RESULT)
            sn = page.evaluate(SNAP_JS)
            T = sn["traces"]
            lines_ = [t for t in T if t["mode"] == "lines" and t["color"]]
            by_style = {(t["color"], t["width"]): t for t in lines_}
            RGB = {"carb": "179,179,230", "lipid": "128,204,204"}
            check(len(T[0]["x"]) == 0 and len(T[1]["x"]) == 0,
                  "SENSITIVITY: on a line map the reaction traces stay EMPTY (the genes' own lines are the reaction arrows; straight "
                  "substrate-to-product lines across the map were the mess in the first version)")
            check(("rgba(%s,1)" % RGB["carb"], 5.5) in by_style and ("rgba(%s,1)" % RGB["lipid"], 5.5) in by_style,
                  f"SENSITIVITY: a HIT line is its pathway CLASS's colour at full strength and thick ({sorted(k for k in by_style if k[1] == 5.5)})")
            check(("rgba(%s,0.65)" % RGB["carb"], 2.5) in by_style and ("rgba(%s,0.32)" % RGB["lipid"], 2) in by_style,
                  "a measured-only line is firmer than an untouched one, both in their class's colour (faint = nobody had it)")
            idx = lambda c, w: by_style[(c, w)]["i"]
            check(idx("rgba(%s,0.32)" % RGB["lipid"], 2) < idx("rgba(%s,0.65)" % RGB["carb"], 2.5) < min(idx("rgba(%s,1)" % RGB["carb"], 5.5), idx("rgba(%s,1)" % RGB["lipid"], 5.5)),
                  "untouched lines are drawn first, then measured, then hits on top")
            hit_carb = by_style[("rgba(%s,1)" % RGB["carb"], 5.5)]
            check(hit_carb["x"][:3] == [0, 200, 400] and hit_carb["y"][:3] == [0, 0, 0],
                  f"a line gene is drawn along its own polyline at its own coordinates, y flipped ({hit_carb['x'][:3]}, {hit_carb['y'][:3]})")
            check(by_style[("rgba(%s,0.32)" % RGB["lipid"], 2)]["y"][:3] == [-300, -300, -300], "...and the untouched one at y = 300")
            decor_t = next(t for t in T if t["mode"] == "lines" and t["color"] == "#E3E8ED")
            check(decor_t["x"] == [0, 400, None] and decor_t["y"] == [-350, -350, None], "KEGG's decoration lines are drawn, faintly")
            chip_ids = {c for t in T if t["custom"] for c in t["custom"]}
            check(chip_ids == {"k1", "k2", "k3", "k5"},
                  f"SENSITIVITY: chips on a HIT or measured line gene and on a hit compound, none on an untouched line ({sorted(chip_ids)})")
            chip_g = next(t for t in T if t["name"] == "Genomics" and t["custom"])
            hand = next(t for t in T if t["mode"] == "markers+text" and t["hover"] and "G12" in " ".join(t["hover"]))
            gx = chip_g["x"][list(chip_g["custom"]).index("k2")] if "k2" in list(chip_g["custom"]) else None
            near = abs(gx - hand["x"][1]) * sn["ppu"] if gx is not None else 1e9
            check(near < 25, f"SENSITIVITY: the chips of a hit line sit at its middle ({near:.1f}px from its handle)")
            check(hand["mcolor"] == "rgba(255,255,255,0)" and "Gmeas" in " ".join(hand["hover"]) and hand["x"][1] == 200,
                  "each line gene has an invisible handle at its middle that carries the hover and the label")
            link_t = next(t for t in T if t["mode"] == "text" and t["text"] is not None and "Glycolysis" in list(t["text"]))
            check(list(link_t["tcolor"]) == ["#B3B3E6"], "a pathway label is printed in its class's colour, as in KEGG's own picture")
            leg = page.evaluate("[...document.querySelectorAll('.multiomics-map-classes input[type=color]')].map(i => [i.getAttribute('data-class'), i.value])")
            check(sorted(leg) == sorted([[CARB, "#b3b3e6"], [LIPID, "#80cccc"]]), f"the class legend has a colour picker per class, set to KEGG's colour ({leg})")
            check("1 hit of 2" in page.inner_text(".multiomics-map-classes"), "each class says how many of its lines were hit")

            print("=== the user recolours a class ===")
            page.locator('.multiomics-map-classes input[data-class="%s"]' % CARB).fill("#ff0000")
            page.wait_for_timeout(1500)
            T2 = page.evaluate(SNAP_JS)["traces"]
            cols2 = {t["color"] for t in T2 if t["mode"] == "lines" and t["color"]}
            check("rgba(255,0,0,1)" in cols2 and "rgba(255,0,0,0.65)" in cols2 and not any("179,179,230" in str(c) for c in cols2) and "rgba(%s,1)" % RGB["lipid"] in cols2,
                  f"SENSITIVITY: choosing red for Carbohydrate recolours exactly that class (hit, measured) and leaves Lipid alone ({sorted(c for c in cols2 if 'rgba' in str(c))})")
            link2 = next(t for t in T2 if t["mode"] == "text" and t["text"] is not None and "Glycolysis" in list(t["text"]))
            check(list(link2["tcolor"]) == ["#ff0000"], "the class's pathway labels follow")
            check(page.evaluate("STATE.multiomics.map.classColors") == {CARB: "#ff0000"}, "the choice is kept for the session")
            check(page.locator("#multiomics-pathway-map").count() == 1, "the chart was replaced, not duplicated")
            page.locator("button.multiomics-map-class-reset").click()
            page.wait_for_timeout(1500)
            T3 = page.evaluate(SNAP_JS)["traces"]
            check(any(t["color"] == "rgba(%s,1)" % RGB["carb"] for t in T3 if t["mode"] == "lines") and page.evaluate("STATE.multiomics.map.classColors") == {},
                  "'KEGG colours' restores the defaults")

            print("=== signalling map: relations, a group ===")
            show(page, SIG_RESULT)
            sn = page.evaluate(SNAP_JS)
            T = sn["traces"]
            rel = [t for t in T if t["mode"] in ("lines+markers", "lines") and t["color"] in ("#1B2A41", "#B5C0CB") and t["name"] == ""]
            kinds = {(t["color"], t["dash"], t["sym"] if isinstance(t["sym"], str) else None): t for t in rel}
            DARK, LIGHT = "#1B2A41", "#B5C0CB"
            check((DARK, "solid", "arrow") in kinds and (LIGHT, "solid", "line-ew") in kinds and (DARK, "solid", None) in kinds
                  and (LIGHT, "dash", "arrow") in kinds and (DARK, "dashdot", None) in kinds,
                  f"SENSITIVITY: activation = arrow, inhibition = bar, binding = plain, indirect = dashed arrow, dissociation = dash-dot; "
                  f"dark where BOTH ends were hit, light otherwise ({sorted((k[0][:4], k[1], k[2] or 'plain') for k in kinds)})")
            n_edges = sum(len([v for v in t["x"] if v is not None]) // 2 for t in rel)
            check(n_edges == 6, f"six relations are drawn; the one pointing at a missing box is skipped ({n_edges})")
            arrow = kinds[(DARK, "solid", "arrow")]
            check(arrow["msize"] == [0, 9, 0, 0, 9, 0], f"only the END of each arrow line carries a head ({arrow['msize']})")
            check(kinds[(DARK, "solid", None)]["msize"] in (0, [0]) or kinds[(DARK, "solid", None)]["mode"] == "lines", "a plain line has no head")
            # the head must stop short of the target's box: AAA(100,100) -> BBB(300,100)
            geo = page.evaluate("""() => { const d = document.getElementById('multiomics-pathway-map'); const xa = d._fullLayout.xaxis, ya = d._fullLayout.yaxis;
                const t = d.data.find(t => t.line && t.line.color === '#1B2A41' && t.marker && t.marker.symbol === 'arrow');
                const b = d._mapPos['b']; const i = t.x.findIndex((x, k) => x !== null && Math.abs(t.y[k] - 100 * -1) < 1e-6 && x > 100 && x < 300);
                return {tx: t.x, ty: t.y, endpx: Math.abs(xa.l2p(t.x[1]) - xa.l2p(b.x)), n: t.x.length}; }""")
            check(6 <= geo["endpx"] <= 40, f"SENSITIVITY: an arrow stops short of its target's box ({geo['endpx']:.1f}px from the centre), so the head is not hidden under it")
            ptm = next(t for t in T if t["mode"] == "text" and t["text"] and "+p" in list(t["text"]))
            check(list(ptm["text"]) == ["+p"], "a phosphorylation is marked '+p' at the middle of its relation")
            groups_t = next(t for t in T if t["mode"] == "lines" and t["color"] == "#A9BDD1")
            xs_ = [v for v in groups_t["x"] if v is not None]
            check(len(groups_t["x"]) == 6 and min(xs_) <= 500 - 6 and max(xs_) >= 600 + 6,
                  f"a group is an outline around ALL its components (x from {min(xs_):.0f} to {max(xs_):.0f}; CCC at 500, EEE at 600)")
            names_ = [t["name"] for t in T if t["name"]]
            check("activation / expression (arrow)" in names_ and "inhibition / repression (bar)" in names_ and "relation with BOTH ends hit (dark)" in names_
                  and "missing interaction (dotted)" not in names_, f"the legend explains the edge styles that are on the map, and only those ({names_})")
            all_hover = " ".join(h for t in T if t["hover"] for h in t["hover"])
            check("CCC+EEE" not in all_hover, "a group is not a box of its own: no hover, no marker, no chips")
            body_ = page.inner_text("body")
            panes = page.locator("button.result-nav-btn").all_inner_texts()
            check(panes == ["Pathway map", "Members", "Interaction links", "Overview"],
                  f"a signalling map has an 'Interaction links' pane and no empty 'Reaction links' one ({panes})")
            page.locator("button.result-nav-btn", has_text="Interaction links").click()
            page.wait_for_timeout(300)
            check("AAA" in page.inner_text("body") and "Genomics; Transcriptomics (up)" in page.inner_text("body"), "the interaction table lists the relation with who hit each end")
            page.locator("button.result-nav-btn", has_text="Overview").click()
            page.wait_for_timeout(300)
            check("Interaction links" in page.inner_text("body"), "the overview counts interaction links")
            page.locator("button.result-nav-btn", has_text="Pathway map").click()
            page.wait_for_timeout(800)

            print("=== a big map is drawn at its own size, in a scroll frame ===")
            show(page, WIDE_RESULT)
            FRAME_JS = """() => { const d = document.getElementById('multiomics-pathway-map'), f = document.querySelector('.multiomics-map-frame');
                return {cw: d.offsetWidth, ch: d.offsetHeight, fw: f.clientWidth, fsw: f.scrollWidth, fh: f.clientHeight, fsh: f.scrollHeight,
                        overflow: getComputedStyle(f).overflow, scrollZoom: d._context && d._context.scrollZoom, sel: document.querySelector('select.multiomics-map-size').value}; }"""
            f = page.evaluate(FRAME_JS)
            check(f["sel"] == "1" and f["cw"] == 2000 + 140 + 20 and f["ch"] == 1200 + 140 + 70,
                  f"SENSITIVITY: a map 2,000 units wide defaults to 100% - drawn {f['cw']} x {f['ch']}px, not squeezed into the page (the unreadable first version)")
            check(f["overflow"] == "auto" and f["fsw"] > f["fw"] + 500 and f["fsh"] > f["fh"] - 1,
                  f"it sits in a scroll frame ({f['fw']}px wide, content {f['fsw']}px)")
            page.evaluate("document.querySelector('.multiomics-map-frame').scrollLeft = 500")
            check(page.evaluate("document.querySelector('.multiomics-map-frame').scrollLeft") == 500, "the frame scrolls sideways")
            check(f["scrollZoom"] is False, "the wheel scrolls the frame (it does not zoom the chart) at a fixed size")
            page.locator("select.multiomics-map-size").select_option("0.5")
            page.wait_for_timeout(1200)
            f2 = page.evaluate(FRAME_JS)
            check(f2["cw"] == round((2000 + 140) * 0.5) + 20, f"the Size menu: 50% is half as wide ({f2['cw']}px)")
            page.locator("select.multiomics-map-size").select_option("2")
            page.wait_for_timeout(1500)
            f3 = page.evaluate(FRAME_JS)
            check(f3["cw"] == (2000 + 140) * 2 + 20, f"200% is twice as wide ({f3['cw']}px)")
            page.locator("select.multiomics-map-size").select_option("fit")
            page.wait_for_timeout(1500)
            f4 = page.evaluate(FRAME_JS)
            check(f4["cw"] <= f4["fw"] and f4["overflow"] == "visible" and f4["scrollZoom"] is True,
                  f"'fit to page' is the page's width again, with no frame and the wheel zooming ({f4['cw']}px in {f4['fw']}px)")
            show(page, LINE_RESULT)
            check(page.locator("select.multiomics-map-size").input_value() == "fit", "a small map defaults to fit-to-page, as before")
            show(page, KEGG_RESULT)
            for size in ("1", "2", "0.5"):
                page.locator("select.multiomics-map-size").select_option(size)
                page.wait_for_timeout(1300)
                mm_ = page.evaluate(CHIP_MIN_JS)
                check(mm_["minDist"] >= MARKER and mm_["three"] >= 1,
                      f"SENSITIVITY: chips on one box stay a marker apart at size {size} ({mm_['minDist']:.1f}px, {mm_['three']} boxes with 3+ chips)")
            show(page, WIDE_RESULT)   # a chart 2,160px wide: its right edge, where the chart's own menu would open, is off-screen
            page.locator("button.multiomics-map-download").click()
            page.wait_for_timeout(400)
            menu = page.evaluate("(() => { const m = document.querySelector('.ora-chart-export-menu'); if (!m) return null; const r = m.getBoundingClientRect(); return {l: r.left, r: r.right, t: r.top, w: innerWidth, items: m.innerText}; })()")
            check(menu is not None and 0 <= menu["l"] and menu["r"] <= menu["w"] and "PNG" in menu["items"],
                  f"'Download map' opens the export menu where you can reach it (the chart's own toolbar is off-screen on a big map): {menu and (round(menu['l']), round(menu['r']))}")
            page.mouse.click(5, 5)

            print("=== the chart's tools, repeated where you can reach them ===")
            show(page, KEGG_RESULT)
            TOOLS_JS = """() => { const d = document.getElementById('multiomics-pathway-map'), fl = d._fullLayout, xa = fl.xaxis, ya = fl.yaxis;
                return {xr: xa.range.slice(), yr: ya.range.slice(), xauto: xa.autorange, dragmode: fl.dragmode, home: d._mapHome,
                        shown: getComputedStyle(document.querySelector('.multiomics-map-tools')).display !== 'none', ppu: Math.abs(xa.l2p(1) - xa.l2p(0))}; }"""
            check(page.evaluate(TOOLS_JS)["shown"] is False, "at 'fit to page' the tools row is hidden (the chart's own toolbar is in view)")
            page.locator("select.multiomics-map-size").select_option("1")
            page.wait_for_timeout(1300)
            t0 = page.evaluate(TOOLS_JS)
            check(t0["shown"] is True and all(page.locator("button.multiomics-map-tool-" + k).count() == 1 for k in ("zoom", "pan", "in", "out", "auto", "reset")),
                  "SENSITIVITY: at a fixed size the chart's toolbar is off-screen, so Zoom tool, Pan tool, Zoom in, Zoom out, Autoscale and Reset axes are offered above it")
            w0 = t0["xr"][1] - t0["xr"][0]
            page.locator("button.multiomics-map-tool-in").click()
            page.wait_for_timeout(800)
            t1 = page.evaluate(TOOLS_JS)
            w1 = t1["xr"][1] - t1["xr"][0]
            check(abs(w1 / w0 - 0.5) < 0.01 and abs((t1["xr"][0] + t1["xr"][1]) / 2 - (t0["xr"][0] + t0["xr"][1]) / 2) < 1e-6 and t1["ppu"] > t0["ppu"] * 1.9,
                  f"SENSITIVITY: Zoom in halves the visible range about its centre ({w1 / w0:.2f} of it, {t1['ppu'] / t0['ppu']:.2f}x the scale)")
            m_in = page.evaluate(CHIP_MIN_JS)
            check(m_in["minDist"] >= MARKER, f"the chips are re-placed for the new zoom and stay a marker apart ({m_in['minDist']:.1f}px)")
            page.locator("button.multiomics-map-tool-out").click()
            page.wait_for_timeout(800)
            t2 = page.evaluate(TOOLS_JS)
            check(abs((t2["xr"][1] - t2["xr"][0]) / w1 - 2) < 0.01, "Zoom out doubles it again")
            page.locator("button.multiomics-map-tool-in").click()
            page.locator("button.multiomics-map-tool-in").click()
            page.wait_for_timeout(800)
            page.locator("button.multiomics-map-tool-reset").click()
            page.wait_for_timeout(800)
            t3 = page.evaluate(TOOLS_JS)
            check(all(abs(a - b) < 1e-6 for a, b in zip(t3["xr"] + t3["yr"], t0["xr"] + t0["yr"])) and (t3["xr"][1] - t3["xr"][0]) > 1.9 * (t1["xr"][1] - t1["xr"][0]),
                  "SENSITIVITY: Reset axes returns exactly to the view the map opened with, however far you zoomed "
                  "(Plotly widens the x range to the canvas's shape, so that - not the raw data box - is the home view)")
            page.locator("button.multiomics-map-tool-in").click()
            page.wait_for_timeout(600)
            page.locator("button.multiomics-map-tool-auto").click()
            page.wait_for_timeout(800)
            t4 = page.evaluate(TOOLS_JS)
            check(t4["xauto"] is True and (t4["xr"][1] - t4["xr"][0]) > w1 * 1.3 and t4["xr"][0] <= 100 and t4["xr"][1] >= 300 and t4["yr"][0] <= -300 and t4["yr"][1] >= -100,
                  f"SENSITIVITY: Autoscale, after zooming in, fits ALL the data again ({t4['xr'][1] - t4['xr'][0]:.0f} units wide, containing every box)")
            page.locator("button.multiomics-map-tool-zoom").click()
            page.wait_for_timeout(500)
            check(page.evaluate(TOOLS_JS)["dragmode"] == "zoom", "the Zoom tool switches a drag to box-zoom")
            page.locator("button.multiomics-map-tool-pan").click()
            page.wait_for_timeout(500)
            check(page.evaluate(TOOLS_JS)["dragmode"] == "pan", "the Pan tool switches it back to pan")
            check(len(errors) == 0, f"no page error from the tools ({errors})")

            print("=== dragging is gone ===")
            show(page, KEGG_RESULT)
            page.evaluate("document.getElementById('multiomics-pathway-map').scrollIntoView({block: 'center'})")
            page.wait_for_timeout(300)
            before = page.evaluate(SNAP_JS)
            a = page.evaluate(PIXELS_JS, "k1")
            page.mouse.move(a["x"], a["y"])
            page.mouse.down()
            page.mouse.move(a["x"] + 40, a["y"] + 30, steps=4)
            page.mouse.up()
            page.wait_for_timeout(500)
            after = page.evaluate(SNAP_JS)
            check(after["pos"] == before["pos"] and True
                  and page.locator("button.multiomics-map-reset").count() == 0,
                  "SENSITIVITY: dragging a box moves nothing, no 'Reset layout' button (the user asked for it to be taken out)")

            print("=== labels ===")
            show(page, BIG_RESULT)
            check(page.locator("select.multiomics-map-labels").input_value() == "hits",
                  "SENSITIVITY: a crowded map (over 300 boxes) labels only the hits by default")
            def labels():
                tt = page.evaluate(SNAP_JS)["traces"]
                return [x for t in tt if t["mode"] == "markers+text" and t["text"] for x in t["text"] if x]
            hits_only = labels()
            check(sorted(hits_only) == sorted(f"Gene{i}" for i in range(320) if i % 50 == 0), f"'hits only' labels the {len(hits_only)} hit boxes and nothing else")
            page.locator("select.multiomics-map-labels").select_option("evidence")
            page.wait_for_timeout(1200)
            ev = labels()
            check(len(ev) == len(hits_only) + len([i for i in range(320) if i % 50 == 1]), f"'hits and measured' adds the measured-only boxes ({len(ev)})")
            page.locator("select.multiomics-map-labels").select_option("all")
            page.wait_for_timeout(1200)
            check(len(labels()) == 320, "'every box' labels all 320")
            page.locator("select.multiomics-map-labels").select_option("none")
            page.wait_for_timeout(1200)
            check(labels() == [], "'none' labels nothing")
            check(page.locator("#multiomics-pathway-map").count() == 1, "re-drawing replaces the chart (no duplicate chart left behind)")
            page.locator("select.multiomics-map-labels").select_option("hits")
            page.wait_for_timeout(1200)
            check(len(errors) == 0, f"no page error through all of this (got {errors})")

            print("=== a chart removed while Plotly is still drawing it ===")
            # A regression check, not a bug fix: the page error first seen in the real-data run came from that test's OWN Plotly.purge
            # of a chart still drawing; the app (which purges the old chart on a label change) raises nothing here, and removing that
            # purge does not change the outcome, so this is NOT claimed to be mutation-sensitive.
            # about the size of KEGG's global map: slow enough to draw that a menu change lands while it is still being drawn
            HUGE_NODES = [dict(node(f"L{i}", "gene", f"L{i}", 10 * (i % 100), 8 * (i // 100), status(genomics=hit(f"L{i}")) if i % 7 == 0 else status(), w=0, h=0),
                               path=[10 * (i % 100) + 3 * k for k in range(1)] and sum(([10 * (i % 100) + 4 * k, 8 * (i // 100) + (k % 3)] for k in range(25)), []))
                          for i in range(1500)] + [node(f"C{i}", "compound", f"C{i}", 7 * (i % 120), 5 * (i // 120), status(), w=8, h=8) for i in range(2500)]
            HUGE = tail("ok", nodes=HUGE_NODES, line_graphics=True)
            page.evaluate("""(res) => { STATE.view = 'annotation_analytics'; STATE.aaPage = 'multiomics_map';
                const mo = freshMultiomicsState(); mo.method = 'mapping'; mo.selected = ['genomics','transcriptomics','proteomics','metabolomics'];
                mo.resolved = {omics: {}}; mo.taxonId = '9606'; mo.taxonName = 'Homo sapiens'; mo.map.result = res; mo.map.database = res.database;
                STATE.multiomics = mo; render(); }""", HUGE)
            n_err = len(errors)
            # no pause at all: switch the labels while the first chart is still being drawn, back and forth
            for _ in range(4):
                page.locator("select.multiomics-map-labels").select_option("all")
                page.locator("select.multiomics-map-labels").select_option("none")
            page.evaluate("""async () => { for (let k = 0; k < 3; k++) { const d = document.createElement('div'); d.style.width = '900px'; document.body.appendChild(d);
                buildPathwayMapPlot(d, STATE.multiomics.map.result, {}); d.remove(); } await new Promise(r => setTimeout(r, 1500)); }""")
            page.wait_for_timeout(6000)
            check(len(errors) == n_err, f"charts torn down while still drawing (huge map, the labels menu switched back and forth with no pause) raise no error ({errors[n_err:][:2]})")
            check(page.locator("#multiomics-pathway-map").count() == 1 and page.evaluate("!!document.getElementById('multiomics-pathway-map').data"),
                  "...and exactly one chart is left, drawn")

            print("=== an overview pathway ===")
            show(page, tail("overview", status="overview", database="Reactome", pathway_id="R-HSA-1643685", pathway_name="Disease",
                            note="Reactome draws \"Disease\" only as links to its 2 sub-pathways.",
                            subpathways={"columns": ["Term ID", "Term Name", "Omics with a hit", "Total hits", "Genomics: hit / annotated", "Joint combined FDR"],
                                         "rows": [{"Term ID": "R-HSA-5668914", "Term Name": "Diseases of metabolism", "Omics with a hit": 2, "Total hits": 9,
                                                   "Genomics: hit / annotated": "9 / 120", "Joint combined FDR": 0.0123},
                                                  {"Term ID": "R-HSA-222", "Term Name": "Quiet", "Omics with a hit": 0, "Total hits": 0,
                                                   "Genomics: hit / annotated": "0 / -", "Joint combined FDR": None}]}))
            body_ = page.inner_text("body")
            check(page.locator("#multiomics-pathway-map").count() == 0 and "only as links to its 2 sub-pathways" in body_,
                  "SENSITIVITY: an overview draws NO chart (the blank map the user hit) and says why")
            check(page.locator("button.multiomics-map-open-sub").count() == 2 and "Diseases of metabolism" in body_ and "9 / 120" in body_ and "0 / -" in body_,
                  "it lists the sub-pathways with their hit counts and a Map button each")
            check(page.locator("button.result-nav-btn").count() == 0, "no empty Members / Links panes on an overview")

        context.close()
    browser.close()

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
