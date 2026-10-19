"""Real-browser test of Multi-omics Cross-omics Pathway Mapping (static/app.js).

Same spirit as playwright_multiomics.py: a UI nobody opened in a browser is not verified. Drives
the whole flow against a live server on port 8009 with real network access (not mocked):
Multi-omics tile screen -> the Cross-omics Pathway Mapping box -> choose all four omics -> one
upload each -> resolve -> database -> list the pathways hit -> search -> open one -> the map, the
member table, the reaction links, the overview -> a Reactome pathway and a sub-pathway without a
diagram of its own -> back to the tile -> re-open (uploads are reused, no re-upload) -> open a
pathway straight from a Joint Pathway Enrichment result.

The things that are easy to get wrong in a UI like this, and are asserted:
  - THE EVIDENCE CHIPS NEVER COVER EACH OTHER. A box hit by several omics carries one chip per
    omics; the check measures the pixel distance between every pair of chips on the same box
    (needs a box with three or more chips, otherwise it would be vacuous) - at first draw AND
    after zooming in and out, because the spacing is in screen pixels while the positions are
    in map units;
  - the map has one trace per omics in that omics' colour, hollow chips for measured-but-not-
    significant, direction triangles for the expression omics, and every hover names the bridge;
  - the tile reuses what is already resolved (no second upload) while Joint still starts clean;
  - a pathway without a diagram says where it was drawn from; downloads are real CSVs;
  - no page error during any of it.
"""
import os
import re

from playwright.sync_api import sync_playwright

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = "http://127.0.0.1:8009"
PASS = FAIL = 0
CHIP_MARKER_PX = 9


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def wait_for(page, expr, seconds, label):
    for _ in range(seconds):
        if page.evaluate(expr):
            return True
        page.wait_for_timeout(1000)
    check(False, f"timed out waiting for: {label}")
    return False


def example(name):
    return os.path.join(_ROOT, "Data3", name)


def pill(page, label):
    return page.locator("label.ora-db-pill", has_text=label).first.locator("input")


def set_omics(page, wanted):
    for label in ["Genomics", "Transcriptomics", "Proteomics", "Metabolomics"]:
        box = pill(page, label)
        if box.is_checked() != (label in wanted):
            box.click()
            page.wait_for_timeout(500)


CHIP_METRICS_JS = """() => {
    const d = document.getElementById('multiomics-pathway-map');
    if (!d || !d._fullLayout || !d.data) return null;
    const xa = d._fullLayout.xaxis, ya = d._fullLayout.yaxis;
    const byNode = {}, names = [], symbols = new Set();
    d.data.forEach(t => {
        if (t.name) names.push(t.name);
        if (!t.customdata || t.mode !== 'markers' || !t.x) return;
        t.customdata.forEach((id, i) => {
            if (typeof id !== 'string' || t.x[i] === null) return;
            (byNode[id] = byNode[id] || []).push([xa.l2p(t.x[i]), ya.l2p(t.y[i])]);
        });
        const sym = t.marker && t.marker.symbol;
        (Array.isArray(sym) ? sym : [sym]).forEach(s => symbols.add(s));
    });
    let minDist = 1e9, multi = 0, three = 0;
    Object.values(byNode).forEach(pts => {
        if (pts.length >= 2) multi++;
        if (pts.length >= 3) three++;
        for (let a = 0; a < pts.length; a++) for (let b = a + 1; b < pts.length; b++)
            minDist = Math.min(minDist, Math.hypot(pts[a][0] - pts[b][0], pts[a][1] - pts[b][1]));
    });
    return {minDist, multi, three, nodes: Object.keys(byNode).length, names, symbols: [...symbols],
            ppu: Math.abs(xa.l2p(1) - xa.l2p(0)), width: d.clientWidth};
}"""


def open_pane(page, title):
    page.locator("button.result-nav-btn", has_text=title).first.click()
    page.wait_for_timeout(600)


def pick_pathway(page, search_text, term_id):
    """Type in the picker's search box and press the Map button of the row with this Term ID."""
    box = page.locator("input.multiomics-map-search")
    box.fill(search_text)
    page.wait_for_timeout(300)
    rows = page.locator(".multiomics-map-picker tbody tr")
    for i in range(rows.count()):
        if term_id in rows.nth(i).inner_text():
            rows.nth(i).locator("button.multiomics-map-open").click()
            return True
    return False


with sync_playwright() as p:
    browser = p.chromium.launch()
    context = browser.new_context(accept_downloads=True, viewport={"width": 1400, "height": 1000})
    page = context.new_page()
    errors = []
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    page.on("console", lambda msg: errors.append("console: " + msg.text) if msg.type == "error" else None)
    page.goto(BASE, wait_until="networkidle")
    page.wait_for_timeout(300)
    check(len(errors) == 0, f"app.js parsed with no top-level error (got {errors})")

    # ----------------------------------------------------------- tile
    print("=== the Cross-omics Pathway Mapping box ===")
    page.evaluate("""() => { STATE.view = 'annotation_analytics'; STATE.aaPage = 'multiomics_home'; render(); }""")
    page.wait_for_timeout(200)
    tiles = page.locator(".aa-genomics-grid button.module-btn")
    check(tiles.nth(1).inner_text().strip() == "Cross-omics Pathway Mapping" and not tiles.nth(1).is_disabled(),
          "the second box is live and shows ONLY the method name (no databases)")
    check(tiles.nth(5).is_disabled() and "coming soon" in tiles.nth(5).inner_text(), "the only unbuilt box (DIABLO) is still greyed out")
    tiles.nth(1).click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "multiomics" and page.evaluate("STATE.multiomics.method") == "mapping",
          "a fresh session opens the upload screen in mapping mode")
    body = page.inner_text("body")
    check("Multi-omics — Cross-omics Pathway Mapping" in body and "Where, inside one pathway" in body,
          "the screen is titled for this method (not Joint Pathway Enrichment)")

    # ----------------------------------------------------------- uploads + resolve
    print("=== four omics: upload and resolve ===")
    set_omics(page, ["Genomics", "Transcriptomics", "Proteomics", "Metabolomics"])
    page.evaluate("""() => { STATE.multiomics.taxonId = '9606'; STATE.multiomics.taxonName = 'Homo sapiens'; render(); }""")
    page.wait_for_timeout(200)
    for i, f in enumerate(["example_multiomics_genes.csv", "example_transcriptomics_de_full.csv",
                           "example_proteomics_de_full.csv", "example_metabolomics_list.csv"]):
        page.locator("input[type=file]").nth(i).set_input_files(example(f))
        page.wait_for_timeout(1800)
    check(page.evaluate("Object.keys(STATE.multiomics.uploads).length") == 4, "all four uploads accepted")
    page.locator("button", has_text="Resolve identifiers").first.click()
    check(wait_for(page, "!!STATE.multiomics.resolved", 1500, "resolve job"), "resolve completed")
    page.wait_for_timeout(300)
    cont = page.locator("button", has_text="Continue to pathway mapping")
    check(cont.count() == 1 and page.locator("button", has_text="Continue to setup").count() == 0,
          "the Continue button leads to the mapping (not to Joint's setup)")
    cont.first.click()
    page.wait_for_timeout(500)
    check(page.evaluate("STATE.aaPage") == "multiomics_map", "navigated to the mapping screen")

    # ----------------------------------------------------------- mapping screen
    print("=== mapping screen: database, cutoffs, picker ===")
    check(page.locator(".ora-config-screen.ora-results-screen").count() == 1, "uses the app's flat results-screen wrapper")
    check(page.locator("select.text-input").first.input_value() == "KEGG", "KEGG is the default database")
    page.wait_for_timeout(1500)
    check("KEGG uses the organism you already chose" in page.inner_text("body"), "KEGG: the shared organism widget is shown")
    check(page.locator("input[type=number]").count() == 4, "four cutoff fields: padj and |log2FC| for each of the two expression omics")
    check(page.locator(".multiomics-map-picker").count() == 0 and page.locator("#multiomics-pathway-map").count() == 0,
          "no picker and no map before the user asks")
    check("computes no p-value" in page.inner_text("body") or "descriptive" in page.inner_text("body"),
          "the screen says up front that this is descriptive")
    page.locator("button", has_text="List pathways hit by my omics").click()
    check(wait_for(page, "!!STATE.multiomics.map.pickers.KEGG", 1500, "KEGG picker"), "the KEGG picker was built")
    page.wait_for_timeout(500)
    picker = page.evaluate("STATE.multiomics.map.pickers.KEGG")
    n_rows = page.locator(".multiomics-map-picker tbody tr").count()
    check(picker["n_pathways"] > 20 and 20 < n_rows <= 150, f"{picker['n_pathways']} pathways listed, {n_rows} rows shown (capped at 150)")
    heads = page.locator(".multiomics-map-picker thead th").all_inner_texts()
    check(heads[-1] == "Map" and "Omics with a hit" in heads and "Joint combined FDR" in heads, f"columns: {heads}")
    check(page.locator(".multiomics-map-picker button.multiomics-map-open").count() == n_rows, "every row has a Map button")
    page.locator("input.multiomics-map-search").fill("citrate cycle")
    page.wait_for_timeout(300)
    n_f = page.locator(".multiomics-map-picker tbody tr").count()
    check(0 < n_f < n_rows and all("itrate cycle" in t for t in page.locator(".multiomics-map-picker tbody tr td:nth-child(2)").all_inner_texts()),
          f"SENSITIVITY: searching narrows the list ({n_rows} -> {n_f}) and every row matches")
    page.locator("input.multiomics-map-search").fill("hsa00020")
    page.wait_for_timeout(300)
    check(page.locator(".multiomics-map-picker tbody tr").count() == 1, "searching by pathway ID finds exactly that pathway")

    # ----------------------------------------------------------- the KEGG map
    print("=== KEGG hsa00020: the map ===")
    page.locator("button.multiomics-map-open").first.click()
    check(wait_for(page, "!!STATE.multiomics.map.result", 600, "mapping job"), "the mapping finished")
    page.wait_for_timeout(1500)
    res = page.evaluate("STATE.multiomics.map.result")
    check(res["status"] == "ok" and res["pathway_id"] == "hsa00020", f"mapped {res['pathway_name']}")
    check("hsa00020" in page.inner_text(".multiomics-map-result, body") and page.locator(".multiomics-map-result").count() == 1,
          "ONE result card, headed with the pathway name and ID")
    for pane in ["Pathway map", "Members", "Reaction links", "Overview"]:
        check(page.locator("button.result-nav-btn", has_text=pane).count() == 1, f"result card has a '{pane}' pane")
    check(page.locator("button.result-nav-btn").first.inner_text() == "Pathway map", "the map is the first pane")
    m = page.evaluate(CHIP_METRICS_JS)
    check(m is not None and m["nodes"] > 10, f"the map is drawn with evidence chips on {m and m['nodes']} boxes")
    check(all(n in m["names"] for n in ["Genomics", "Transcriptomics", "Proteomics", "Metabolomics"]),
          f"one trace per omics, named (legend: {[n for n in m['names'] if n]})")
    check("circle" in m["symbols"] and "circle-open" in m["symbols"],
          f"chips: circles for list hits and hollow circles for measured-only ({sorted(set(map(str, m['symbols'])))}) - the synthetic "
          f"expression tables plant no TCA signal, so no triangles are expected on THIS pathway (checked on another below)")
    check(m["three"] >= 1,
          f"the overlap case exists: {m['three']} boxes carry three or more chips (without this the next check would be vacuous)")
    check(m["minDist"] >= CHIP_MARKER_PX,
          f"SENSITIVITY: chips on the same box are at least a marker-width apart at first draw (closest pair {m['minDist']:.1f}px; "
          f"a marker is {CHIP_MARKER_PX}px - stacked, this would be 0)")
    # zoom in, then out: spacing is in screen pixels, so it must hold at every zoom level
    page.evaluate("""() => { const d = document.getElementById('multiomics-pathway-map');
        const r = d._fullLayout.xaxis.range; const mid = (r[0] + r[1]) / 2, w = (r[1] - r[0]) / 5;
        Plotly.relayout(d, {'xaxis.range': [mid - w, mid + w]}); }""")
    page.wait_for_timeout(1200)
    mz = page.evaluate(CHIP_METRICS_JS)
    check(mz["ppu"] > m["ppu"] * 2 and mz["minDist"] >= CHIP_MARKER_PX,
          f"SENSITIVITY: zoomed in ({m['ppu']:.2f} -> {mz['ppu']:.2f} px per map unit) chips stay clear ({mz['minDist']:.1f}px)")
    page.evaluate("""() => { const d = document.getElementById('multiomics-pathway-map');
        Plotly.relayout(d, {'xaxis.autorange': true, 'yaxis.autorange': true}); }""")
    page.wait_for_timeout(1200)
    mo = page.evaluate(CHIP_METRICS_JS)
    check(mo["minDist"] >= CHIP_MARKER_PX, f"...and after zooming back out ({mo['minDist']:.1f}px)")
    hover = page.evaluate("""() => { const d = document.getElementById('multiomics-pathway-map');
        const t = d.data.find(t => t.name === 'Genomics' && t.hovertext); return t ? t.hovertext[0] : null; }""")
    check(hover and "Genomics: hit" in hover and "[" in hover and "<b>" in hover,
          f"a chip's hover shows the box, the state and the identifier bridge ({(hover or '')[:90]!r})")
    node_hover = page.evaluate("""() => { const d = document.getElementById('multiomics-pathway-map');
        const t = d.data.find(t => t.mode === 'markers+text' && t.hovertext.some(h => h.includes('Genomics: hit')));
        return t ? t.hovertext.find(h => h.includes('Genomics: hit')) : null; }""")
    check(node_hover is not None and "<b>" in node_hover, "a box's hover lists every omics that hit it")
    check("filled = hit, hollow = measured but not significant" in page.inner_text("body"), "the pane explains the chips")

    # ----------------------------------------------------------- members / links / overview
    print("=== members, reaction links, overview ===")
    open_pane(page, "Members")
    check(page.locator(".multiomics-map-result table:visible tbody tr").count() >= 20, "the member table has rows")
    heads = page.locator(".multiomics-map-result table:visible thead th").all_inner_texts()
    check(heads == ["Pathway node", "Box ID", "Node type", "Member", "Omics", "Status", "Direction", "log2FC", "padj", "Matched by"], f"member columns: {heads}")
    texts = page.locator(".multiomics-map-result table:visible tbody").first.inner_text()
    check("measured, not significant" in texts and "UniProt -> KEGG gene" in texts,
          "SENSITIVITY: measured-not-significant members are listed beside hits, with how each was matched")
    with page.expect_download() as dl:
        page.locator("button", has_text="Download members (CSV)").click()
    path = dl.value.path()
    csv = open(path, encoding="utf-8").read().splitlines()
    check(csv[0].startswith("Pathway node,Box ID,Node type,Member,Omics") and len(csv) - 1 >= 20, f"the members CSV downloads ({len(csv) - 1} rows)")
    open_pane(page, "Reaction links")
    link_text = page.locator(".multiomics-map-result table:visible tbody").first.inner_text()
    check("itrate" in link_text and "CS" in link_text, "BIOLOGY: citrate synthase (CS) is linked to citrate")
    with page.expect_download() as dl2:
        page.locator("button", has_text="Download reaction links (CSV)").click()
    check(open(dl2.value.path(), encoding="utf-8").read().startswith("Reaction,Enzyme (hit)"), "the links CSV downloads")
    open_pane(page, "Overview")
    body = page.inner_text("body")
    check("computes no p-value" in body and "Boxes hit by any omics" in body and "Cross-omics reaction links" in body,
          "overview: the descriptive caveat and the metric cards")
    check(page.locator(".multiomics-map-result table").count() >= 1 and "Boxes measured only" in body, "overview: per-omics table with the measured-only count")
    open_pane(page, "Pathway map")

    # ----------------------------------------------------------- a pathway the expression omics DO hit
    print("=== a pathway the expression tables hit: direction triangles ===")
    best = page.evaluate("""() => { const rows = STATE.multiomics.map.pickers.KEGG.table.rows;
        let best = null, bk = -1;
        rows.forEach(r => { const k = parseInt(String(r['Transcriptomics: hit / annotated']).split(' / ')[0]) +
                                      parseInt(String(r['Proteomics: hit / annotated']).split(' / ')[0]);
                            if (k > bk && r['Total hits'] < 400) { bk = k; best = r['Term ID']; } });
        return [best, bk]; }""")
    check(best[1] >= 5, f"the picker has a pathway with {best[1]} expression-omics hits ({best[0]})")
    check(pick_pathway(page, best[0], best[0]) and wait_for(
        page, f"STATE.multiomics.map.result && STATE.multiomics.map.result.pathway_id === '{best[0]}'", 600, "DE pathway"),
        "mapped it")
    page.wait_for_timeout(2000)
    de = page.evaluate(CHIP_METRICS_JS)
    check("triangle-up" in de["symbols"] or "triangle-down" in de["symbols"],
          f"expression-omics hits are triangles ({sorted(set(map(str, de['symbols'])))})")
    consistent = page.evaluate("""() => { const d = document.getElementById('multiomics-pathway-map'); let n = 0, bad = 0;
        d.data.forEach(t => { if (!['Transcriptomics', 'Proteomics'].includes(t.name) || !t.hovertext) return;
            const sym = t.marker.symbol; t.hovertext.forEach((h, i) => {
                const up = h.includes('hit (up)'), down = h.includes('hit (down)'), mixed = h.includes('hit (mixed)');
                const s = Array.isArray(sym) ? sym[i] : sym; n++;
                if ((up && s !== 'triangle-up') || (down && s !== 'triangle-down') || (mixed && s !== 'diamond')) bad++; }); });
        return [n, bad]; }""")
    check(consistent[0] >= 3 and consistent[1] == 0,
          f"SENSITIVITY: every expression chip's shape agrees with the direction in its own hover ({consistent[0]} chips, {consistent[1]} disagree)")
    hov = page.evaluate("""() => { const d = document.getElementById('multiomics-pathway-map');
        const t = d.data.find(t => t.name === 'Transcriptomics' && t.hovertext); return t ? t.hovertext[0] : null; }""")
    check(hov and "log2FC" in hov and "padj" in hov and "[" in hov, f"an expression chip's hover shows log2FC, padj and the bridge ({(hov or '')[:100]!r})")
    check(de["minDist"] >= CHIP_MARKER_PX, f"chips clear on this denser map too ({de['minDist']:.1f}px, {de['nodes']} boxes)")

    # ----------------------------------------------------------- KEGG overview map
    print("=== KEGG's global map hsa01100: genes are lines ===")
    check(pick_pathway(page, "hsa01100", "hsa01100") and wait_for(
        page, "STATE.multiomics.map.result && STATE.multiomics.map.result.pathway_id === 'hsa01100'", 900, "hsa01100"),
        "mapped Metabolic pathways (3,562 compounds, 1,342 genes)")
    page.wait_for_timeout(2500)
    g = page.evaluate("""() => { const d = document.getElementById('multiomics-pathway-map'); if (!d || !d.data) return null;
        const xa = d._fullLayout.xaxis, ya = d._fullLayout.yaxis, T = d.data;
        const handle = T.find(t => t.mode === 'markers+text' && t.marker && t.marker.color === 'rgba(255,255,255,0)');
        let minx = Infinity, miny = Infinity; handle.x.forEach((x, i) => { minx = Math.min(minx, x); miny = Math.min(miny, handle.y[i]); });
        let near = 0; handle.x.forEach((x, i) => { if (Math.hypot(xa.l2p(x) - xa.l2p(minx), ya.l2p(handle.y[i]) - ya.l2p(miny)) < 4) near++; });
        const polylines = (t) => t.x.filter(v => v === null).length;
        const lineTr = T.filter(t => t.mode === 'lines' && t.line && /^rgba\(\d+,\d+,\d+,/.test(t.line.color));
        const byWidth = (w) => lineTr.filter(t => t.line.width === w).reduce((a, t) => a + polylines(t), 0);
        const faint = byWidth(2) + byWidth(2.5), thick = byWidth(5.5);
        const labels = T.filter(t => t.mode === 'markers+text' && t.text).reduce((a, t) => a + t.text.filter(x => x).length, 0);
        const chipIds = new Set(); T.forEach(t => { if (t.customdata && t.mode === 'markers') t.customdata.forEach(c => chipIds.add(c)); });
        const res = STATE.multiomics.map.result, genes = res.nodes.filter(n => n.kind === 'gene');
        const hitGenes = genes.filter(n => Object.values(n.status).some(v => v.state === 'hit')).length;
        const hitLineIds = genes.filter(n => Object.values(n.status).some(v => v.state === 'hit')).map(n => n.id);
        return {handles: handle.x.length, near, faint, thick, genes: genes.length, hitGenes, chipsOnHits: hitLineIds.filter(id => chipIds.has(id)).length,
                labels, nodes: res.nodes.length, labelSel: document.querySelector('select.multiomics-map-labels').value, reactionPts: T[0].x.length,
                decor: T.some(t => t.mode === 'lines' && t.line && t.line.color === '#E3E8ED'),
                sizeSel: document.querySelector('select.multiomics-map-size').value, chartW: d.offsetWidth, chartH: d.offsetHeight,
                frameW: document.querySelector('.multiomics-map-frame').clientWidth, classes: res.classes.length,
                swatches: document.querySelectorAll('.multiomics-map-classes input[type=color]').length}; }""")
    check(g is not None and g["handles"] == 1342 and g["genes"] == 1342, f"every one of the {g and g['genes']} genes is on the chart (one handle each)")
    check(g["near"] <= 3, f"SENSITIVITY: the genes are NOT piled in a corner ({g['near']} within 4px of the top-left-most gene; before: all 1,342)")
    check(g["thick"] == g["hitGenes"] and g["hitGenes"] > 20 and g["faint"] == g["genes"] - g["hitGenes"],
          f"SENSITIVITY: the {g['hitGenes']} hit genes are thick full-strength lines, the other {g['faint']} faint or firm - none lost or doubled")
    check(g["chipsOnHits"] == g["hitGenes"], f"every hit line carries its omics chips ({g['chipsOnHits']} of {g['hitGenes']})")
    check(g["sizeSel"] == "1" and g["chartW"] >= 4000 and g["frameW"] < 1400 and g["chartH"] >= 2500,
          f"SENSITIVITY: the map is drawn at its own size ({g['chartW']} x {g['chartH']}px) in a frame {g['frameW']}px wide that you scroll - not squeezed into the page")
    check(g["classes"] >= 8 and g["swatches"] == g["classes"], f"the class legend has {g['swatches']} colour pickers for the {g['classes']} pathway classes")
    check(g["reactionPts"] == 0 and g["decor"], "straight substrate-to-product lines are not drawn on a line map; KEGG's decoration lines are")
    hit_boxes = page.evaluate("STATE.multiomics.map.result.nodes.filter(n => n.kind !== 'link' && n.kind !== 'group' && Object.values(n.status).some(v => v.state === 'hit')).length")
    check(g["labelSel"] == "hits" and g["labels"] == hit_boxes and hit_boxes < g["nodes"] / 10,
          f"SENSITIVITY: with {g['nodes']} boxes only the hits are labelled by default ({g['labels']} labels = the {hit_boxes} hit boxes, not {g['nodes']})")
    t_draw = page.evaluate("""async () => { const t0 = performance.now(); const d = document.createElement('div'); d.style.width = '1000px';
        document.body.appendChild(d); buildPathwayMapPlot(d, STATE.multiomics.map.result, {labels: 'hits'});
        // _mapReady is set by the chart's onReady, i.e. once Plotly has FINISHED drawing (_fullLayout exists long before that)
        while (!d._mapReady) await new Promise(r => setTimeout(r, 20));
        const dt = performance.now() - t0; d.remove(); return dt; }""")
    check(t_draw < 8000, f"the biggest map (about 5,000 boxes) is drawn in {t_draw / 1000:.1f} s once its data is here (timed to the chart's own draw-complete flag; the first request for it also downloads KEGG's file)")
    check(not any(e for e in errors), f"no error drawing it ({errors})")
    page.locator("select.multiomics-map-labels").select_option("all")
    page.wait_for_timeout(3000)
    n_all = page.evaluate("document.getElementById('multiomics-pathway-map').data.filter(t => t.mode === 'markers+text' && t.text).reduce((a, t) => a + t.text.filter(x => x).length, 0)")
    check(n_all > g["labels"] * 5, f"'every box' labels far more ({n_all} labels)")
    page.locator("select.multiomics-map-labels").select_option("hits")
    page.wait_for_timeout(3000)
    print("=== the chart's tools on the real global map ===")
    TOOLS = """() => { const d = document.getElementById('multiomics-pathway-map'), fl = d._fullLayout; return {w: fl.xaxis.range[1] - fl.xaxis.range[0], ppu: Math.abs(fl.xaxis.l2p(1) - fl.xaxis.l2p(0))}; }"""
    t_a = page.evaluate(TOOLS)
    page.locator("button.multiomics-map-tool-in").click()
    page.wait_for_timeout(2500)
    t_b = page.evaluate(TOOLS)
    check(page.locator(".multiomics-map-tools").is_visible() and abs(t_b["w"] / t_a["w"] - 0.5) < 0.01 and t_b["ppu"] > t_a["ppu"] * 1.9,
          f"SENSITIVITY: the zoom / pan / autoscale / reset tools are on the page for the 5,000px map and Zoom in works on it ({t_b['w'] / t_a['w']:.2f} of the range)")
    page.locator("button.multiomics-map-tool-reset").click()
    page.wait_for_timeout(2500)
    check(abs(page.evaluate(TOOLS)["w"] / t_a["w"] - 1) < 1e-6, "Reset axes returns to the whole map")

    print("=== the user recolours a pathway class ===")
    cls = page.evaluate("STATE.multiomics.map.result.classes")
    lip = next(c for c in cls if c["name"] == "Lipid metabolism")
    page.locator('.multiomics-map-classes input[data-class="Lipid metabolism"]').fill("#ff0000")
    page.wait_for_timeout(5000)
    red = page.evaluate("""() => { const d = document.getElementById('multiomics-pathway-map');
        return d.data.filter(t => t.mode === 'lines' && t.line && /^rgba\\(255,0,0,/.test(t.line.color)).reduce((a, t) => a + t.x.filter(v => v === null).length, 0); }""")
    check(red == lip["n_lines"], f"SENSITIVITY: choosing red for Lipid metabolism recolours exactly its {lip['n_lines']} lines ({red} are red)")
    page.locator("button.multiomics-map-class-reset").click()
    page.wait_for_timeout(5000)

    print("=== labels need room: measured overlap at the page's width and at 100% ===")
    OVERLAP_JS = """() => { const d = document.getElementById('multiomics-pathway-map'); const xa = d._fullLayout.xaxis, ya = d._fullLayout.yaxis; const boxes = [];
        d.data.forEach(t => { if (t.mode !== 'markers+text' || !t.text) return; t.text.forEach((s, i) => { if (!s) return;
            const x = xa._offset + xa.l2p(t.x[i]), y = ya._offset + ya.l2p(t.y[i]), w = s.length * 5.6; boxes.push([x - w / 2, y - 15, x + w / 2, y - 4]); }); });
        let over = 0; for (let i = 0; i < boxes.length; i++) for (let j = 0; j < boxes.length; j++)
            if (i !== j && boxes[i][0] < boxes[j][2] && boxes[j][0] < boxes[i][2] && boxes[i][1] < boxes[j][3] && boxes[j][1] < boxes[i][3]) { over++; break; }
        return {n: boxes.length, over}; }"""
    page.locator("select.multiomics-map-size").select_option("fit")
    page.wait_for_timeout(6000)
    ov_fit = page.evaluate(OVERLAP_JS)
    page.locator("select.multiomics-map-size").select_option("1")
    page.wait_for_timeout(6000)
    ov_100 = page.evaluate(OVERLAP_JS)
    f_fit, f_100 = ov_fit["over"] / max(ov_fit["n"], 1), ov_100["over"] / max(ov_100["n"], 1)
    print(f"      measured: {ov_fit['over']} of {ov_fit['n']} hit labels overlap another when squeezed to the page ({f_fit:.0%}); {ov_100['over']} of {ov_100['n']} at 100% ({f_100:.0%})")
    check(ov_100["n"] == ov_fit["n"] and f_100 < f_fit and f_100 < 0.35,
          f"SENSITIVITY: drawn at its own size far fewer labels sit on top of each other ({f_100:.0%} at 100% against {f_fit:.0%} squeezed to the page)")

    # ----------------------------------------------------------- KEGG signalling map
    print("=== KEGG's signalling map hsa04210 (apoptosis): relations and groups ===")
    check(pick_pathway(page, "hsa04210", "hsa04210") and wait_for(
        page, "STATE.multiomics.map.result && STATE.multiomics.map.result.pathway_id === 'hsa04210'", 600, "hsa04210"), "mapped Apoptosis")
    page.wait_for_timeout(2500)
    sg = page.evaluate("""() => { const d = document.getElementById('multiomics-pathway-map'); const T = d.data, res = STATE.multiomics.map.result;
        const rel = T.filter(t => t.mode === 'lines+markers' || (t.mode === 'lines' && t.line && ['#1B2A41', '#B5C0CB'].includes(t.line.color)));
        const edges = rel.reduce((a, t) => a + t.x.filter(v => v !== null).length / 2, 0);
        const dark = T.filter(t => t.line && t.line.color === '#1B2A41').reduce((a, t) => a + t.x.filter(v => v !== null).length / 2, 0);
        const groups = T.find(t => t.line && t.line.color === '#A9BDD1');
        const hit = (n) => Object.values(n.status).some(v => v.state === 'hit'); const by = {}; res.nodes.forEach(n => by[n.id] = n);
        const bothHit = res.relations.filter(r => by[r.a] && by[r.b] && hit(by[r.a]) && hit(by[r.b])).length;
        return {edges, dark, bothHit, nrel: res.relations.length, groupRects: groups ? groups.x.filter(v => v === null).length : 0,
                ngroups: res.nodes.filter(n => n.kind === 'group').length, legend: T.map(t => t.name).filter(Boolean),
                panes: [...document.querySelectorAll('button.result-nav-btn')].map(b => b.innerText), nint: res.summary.n_interactions}; }""")
    check(sg["nrel"] == 135 and sg["edges"] >= 120, f"the map's own {sg['nrel']} relations are drawn as lines ({sg['edges']:.0f} edges; very short ones are left out)")
    check(sg["dark"] == sg["bothHit"] or abs(sg["dark"] - sg["bothHit"]) <= 8 and sg["bothHit"] > 0,
          f"SENSITIVITY: the dark edges are the relations with both ends hit ({sg['dark']:.0f} drawn, {sg['bothHit']} by recount)")
    check(sg["ngroups"] == 4 and sg["groupRects"] == 4, f"the {sg['ngroups']} groups are outlined")
    check(any("arrow" in n for n in sg["legend"]) and any("bar" in n for n in sg["legend"]), f"the legend explains arrows and bars ({sg['legend'][-6:]})")
    check(sg["panes"] == ["Pathway map", "Members", "Interaction links", "Overview"], f"panes: {sg['panes']}")
    page.locator("button.result-nav-btn", has_text="Interaction links").click()
    page.wait_for_timeout(400)
    check(page.locator(".multiomics-map-result table:visible tbody tr").count() == sg["nint"] and sg["nint"] > 0,
          f"the interaction table has one row per interaction link ({sg['nint']})")
    with page.expect_download() as dli:
        page.locator("button", has_text="Download interaction links (CSV)").click()
    check(open(dli.value.path(), encoding="utf-8").read().startswith("From,To,Relation,From hit by,To hit by"), "the interaction links download as CSV")
    page.locator("button.result-nav-btn", has_text="Pathway map").click()
    page.wait_for_timeout(1500)
    chips_sg = page.evaluate("""() => { const d = document.getElementById('multiomics-pathway-map'); const res = STATE.multiomics.map.result;
        const ids = new Set(); d.data.forEach(t => { if (t.customdata && t.mode === 'markers') t.customdata.forEach(c => ids.add(c)); });
        const hit = res.nodes.filter(n => n.kind !== 'link' && n.kind !== 'group' && Object.values(n.status).some(v => v.state === 'hit'));
        return {hit: hit.length, withChips: hit.filter(n => ids.has(n.id)).length, sel: document.querySelector('select.multiomics-map-size').value}; }""")
    check(chips_sg["hit"] > 20 and chips_sg["withChips"] == chips_sg["hit"], f"every hit box of the signalling map carries its chips ({chips_sg['withChips']} of {chips_sg['hit']})")

    # ----------------------------------------------------------- Reactome
    print("=== Reactome: a diagram, and a sub-pathway drawn from its ancestor ===")
    page.locator("select.text-input").first.select_option("Reactome")
    page.wait_for_timeout(400)
    check("Reactome draws the pathway from its own diagram" in page.inner_text("body"), "Reactome: its own explanation replaces the KEGG organism widget")
    page.locator("button", has_text="List pathways hit by my omics").click()
    check(wait_for(page, "!!STATE.multiomics.map.pickers.Reactome", 1500, "Reactome picker"), "the Reactome picker was built")
    page.wait_for_timeout(500)
    check(pick_pathway(page, "R-HSA-71403", "R-HSA-71403"), "found the citric acid cycle in the picker")
    check(wait_for(page, "!!STATE.multiomics.map.result", 600, "Reactome mapping"), "Reactome pathway mapped")
    page.wait_for_timeout(1500)
    rres = page.evaluate("STATE.multiomics.map.result")
    rm = page.evaluate(CHIP_METRICS_JS)
    check(rres["database"] == "Reactome" and rm is not None and rm["nodes"] >= 5, f"the Reactome map is drawn ({rm and rm['nodes']} boxes with chips)")
    check(rm["minDist"] >= CHIP_MARKER_PX, f"chips clear on the Reactome diagram too ({rm['minDist']:.1f}px)")
    check(pick_pathway(page, "R-HSA-70171", "R-HSA-70171") and wait_for(page, "STATE.multiomics.map.result && STATE.multiomics.map.result.pathway_id === 'R-HSA-70171'", 600, "glycolysis"),
          "Glycolysis (no diagram of its own) mapped")
    page.wait_for_timeout(1500)
    open_pane(page, "Overview")
    check("has no diagram of its own" in page.inner_text("body") and "Glucose metabolism" in page.inner_text("body"),
          "the overview says where a diagram-less pathway was drawn from")
    open_pane(page, "Pathway map")


    print("=== a Reactome overview pathway (Disease) ===")
    check(pick_pathway(page, "R-HSA-1643685", "R-HSA-1643685") and wait_for(
        page, "STATE.multiomics.map.result && STATE.multiomics.map.result.pathway_id === 'R-HSA-1643685'", 600, "Disease"), "opened Disease")
    page.wait_for_timeout(1500)
    check(page.evaluate("STATE.multiomics.map.result.status") == "overview" and page.locator("#multiomics-pathway-map").count() == 0,
          "SENSITIVITY: Disease shows NO chart - the blank map the user hit")
    body = page.inner_text("body")
    check("only as links to its 12 sub-pathways" in body and "Diseases of metabolism" in body,
          "it says why and lists the sub-pathways")
    depth = 0
    while page.locator("#multiomics-pathway-map").count() == 0 and depth < 5 and page.locator("button.multiomics-map-open-sub").count():
        page.locator("button.multiomics-map-open-sub").first.click()
        wait_for(page, "!!STATE.multiomics.map.result", 600, "sub-pathway")
        page.wait_for_timeout(2500)
        depth += 1
    check(page.locator("#multiomics-pathway-map").count() == 1 and page.evaluate("STATE.multiomics.map.result.status") == "ok",
          f"SENSITIVITY: pressing Map on the sub-pathways reaches a real diagram after {depth} step(s): "
          f"{page.evaluate('STATE.multiomics.map.result.pathway_name')}")
    check(page.evaluate("STATE.multiomics.map.result.nodes.length") > 5, "...with boxes on it")

    # ----------------------------------------------------------- back / reuse
    print("=== navigation and reuse of the uploads ===")
    check(page.locator("button", has_text="Back to upload").count() == 1 and page.locator("button", has_text="Back to Multi-omics").count() == 1,
          "the mapping screen has both ways back")
    page.locator("button", has_text="Back to Multi-omics").click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "multiomics_home", "Back to Multi-omics returns to the tile screen")
    page.locator(".aa-genomics-grid button.module-btn").nth(1).click()
    page.wait_for_timeout(500)
    check(page.evaluate("STATE.aaPage") == "multiomics_map" and page.evaluate("!!STATE.multiomics.resolved"),
          "SENSITIVITY: re-opening the box goes straight to the picker - the resolved uploads are reused, not asked for again")

    # ----------------------------------------------------------- from a joint result
    print("=== open a pathway from a Joint Pathway Enrichment result ===")
    page.locator("button", has_text="Back to Multi-omics").click()
    page.wait_for_timeout(300)
    page.locator(".aa-genomics-grid button.module-btn").nth(0).click()
    page.wait_for_timeout(400)
    check(page.evaluate("STATE.aaPage") == "multiomics" and page.evaluate("STATE.multiomics.method") == "joint"
          and page.evaluate("!STATE.multiomics.resolved"), "the Joint box still starts from a clean setup")
    set_omics(page, ["Genomics", "Metabolomics"])
    page.evaluate("""() => { STATE.multiomics.taxonId = '9606'; STATE.multiomics.taxonName = 'Homo sapiens'; render(); }""")
    page.wait_for_timeout(200)
    page.locator("input[type=file]").nth(0).set_input_files(example("example_multiomics_genes.csv"))
    page.wait_for_timeout(1500)
    page.locator("input[type=file]").nth(1).set_input_files(example("example_metabolomics_list.csv"))
    page.wait_for_timeout(1500)
    page.locator("button", has_text="Resolve identifiers").first.click()
    check(wait_for(page, "!!STATE.multiomics.resolved", 900, "joint resolve"), "Genomics + Metabolomics resolved")
    page.locator("button", has_text="Continue to setup").first.click()
    page.wait_for_timeout(500)
    pill(page, "Reactome").click()
    page.wait_for_timeout(300)
    page.locator("button.ora-run-btn").click()
    check(wait_for(page, "!!STATE.multiomics.run.results", 900, "joint run"), "a KEGG joint run completed")
    page.wait_for_timeout(600)
    open_pane(page, "Combined table")
    check(page.locator("select.multiomics-map-from-joint").count() == 1 and page.locator("button.multiomics-map-from-joint-btn").count() == 1,
          "the combined table offers 'Map this pathway'")
    first_value = page.locator("select.multiomics-map-from-joint").input_value()
    page.locator("button.multiomics-map-from-joint-btn").click()
    check(wait_for(page, "STATE.aaPage === 'multiomics_map' && !!STATE.multiomics.map.result", 600, "map from joint"),
          "clicking it opens the mapping screen with that pathway already mapped")
    page.wait_for_timeout(1500)
    check(page.evaluate("STATE.multiomics.map.result.pathway_id") == first_value and page.evaluate("STATE.multiomics.method") == "mapping",
          f"it is the pathway that was chosen ({first_value})")
    check(page.evaluate("STATE.multiomics.map.result.omics") == ["genomics", "metabolomics"], "mapped with the omics of that joint run")
    check(page.evaluate(CHIP_METRICS_JS) is not None, "the map is drawn")
    page.locator("button", has_text="List pathways hit by my omics").click()
    check(wait_for(page, "!!STATE.multiomics.map.pickers.KEGG", 900, "picker after joint"), "the picker works afterwards")
    page.wait_for_timeout(400)
    check("Joint Pathway Enrichment FDR" in page.inner_text("body"), "the picker shows the Joint FDR from the earlier run")

    print("=== console ===")
    check(len(errors) == 0, f"no page error or console error during the whole flow (got {errors})")
    browser.close()

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
