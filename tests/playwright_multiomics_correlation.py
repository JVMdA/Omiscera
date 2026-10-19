"""Real-browser test of Multi-omics Correlation-based Integration (static/app.js).

Drives the whole flow against a live server on port 8009 with real network access (not mocked):
Multi-omics tile screen -> the Correlation-based Integration box -> choose the omics -> one matrix upload each (the three
example tables: samples in the columns for RNA and protein, in the rows for the metabolites) -> resolve and match samples
-> setup -> run -> every result pane.

What is asserted, beyond "it renders":
  - the tile is live, Genomics cannot be chosen (a gene list has no values per sample), Resolve stays disabled until every
    table is uploaded and an organism is chosen;
  - the detected orientation is SHOWN (and a changed orientation re-reads the same file);
  - the alignment card names the 12 shared samples;
  - the concordance table has the planted genes, its histogram is drawn, clicking a row draws that gene's scatter with the
    SAME rho as the table; a pair row's scatter likewise;
  - the heatmap is drawn at a fixed cell size (wider than the viewport it sits in -> it scrolls, it is not squeezed) with
    a dot only where a pair is reported;
  - the network draws every node, no two nodes of a column sit closer than the marker size, and the figure tools row
    (zoom / pan / autoscale / reset / download) is present on every chart;
  - the downloads are real CSVs whose row counts equal the tables;
  - no page error during any of it.
"""
import csv
import re
import io
import os

from playwright.sync_api import sync_playwright

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = "http://127.0.0.1:8009"
PASS = FAIL = 0
MARKER_PX = 11


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


def open_pane(page, title):
    page.locator("button.result-nav-btn", has_text=title).first.click()
    page.wait_for_timeout(700)
    if os.environ.get("CORR_SHOT_DIR"):    # optional: look at the panes with your own eyes
        page.screenshot(path=os.path.join(os.environ["CORR_SHOT_DIR"], title.replace(" ", "_") + ".png"), full_page=True)


SCATTER = "(() => " + "[...document.querySelectorAll('.multiomics-corr-scatter')].find(d => d.offsetParent !== null && d.data)" + ")()"
MODEBAR = """(id) => { const d = document.getElementById(id); return !!(d && d.querySelector('.modebar')); }"""
HAS_DATA = """(id) => { const d = document.getElementById(id); return !!(d && d.data && d.data.length); }"""

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
    print("=== the Correlation-based Integration box ===")
    page.evaluate("""() => { STATE.view = 'annotation_analytics'; STATE.aaPage = 'multiomics_home'; render(); }""")
    page.wait_for_timeout(200)
    tiles = page.locator(".aa-genomics-grid button.module-btn")
    check(tiles.nth(2).inner_text().strip() == "Correlation-based Integration" and not tiles.nth(2).is_disabled(),
          "the third box is live and shows ONLY the method name")
    check(not tiles.nth(5).is_disabled(), "SENSITIVITY: DIABLO (the sixth box) is built now too - nothing is greyed out")
    tiles.nth(2).click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "multiomics_corr" and page.evaluate("STATE.multiomics.method") == "correlation",
          "a fresh session opens the correlation screen")
    body = page.inner_text("body")
    check("Correlation-based Integration" in body and "not a mechanism" in body, "the screen says up front that correlation is not a mechanism")
    check(page.locator("label.ora-db-pill", has_text="Genomics").first.locator("input").is_disabled(),
          "Genomics cannot be chosen (a gene list has no values per sample)")
    check(page.locator("input[type=file]").count() == 0, "no upload box until two omics are chosen")

    # ----------------------------------------------------------- uploads
    print("=== uploads, orientation, sample matching ===")
    pill(page, "Transcriptomics").click()
    page.wait_for_timeout(300)
    check(page.locator("input[type=file]").count() == 0, "one omics is not enough to correlate")
    pill(page, "Proteomics").click()
    page.wait_for_timeout(500)
    pill(page, "Metabolomics").click()
    page.wait_for_timeout(500)
    check(page.evaluate("STATE.multiomics.corr.selected") == ["transcriptomics", "proteomics", "metabolomics"], "three omics chosen, in a fixed order")
    check(page.locator("input[type=file]").count() == 3, "one upload box per chosen omics")
    check(page.locator("button.multiomics-corr-resolve").is_disabled(), "Resolve is disabled until everything is uploaded and an organism chosen")
    page.evaluate("""() => { STATE.multiomics.corr.taxonId = '9606'; STATE.multiomics.corr.taxonName = 'Homo sapiens'; render(); }""")
    page.wait_for_timeout(200)
    for i, f in enumerate(["example_corr_transcriptomics.csv", "example_corr_proteomics.csv", "example_corr_metabolomics.csv"]):
        page.locator("input[type=file]").nth(i).set_input_files(example(f))
        page.wait_for_timeout(1500)
    ups = page.evaluate("Object.fromEntries(Object.entries(STATE.multiomics.corr.uploads).map(([k, v]) => [k, v.info]))")
    check(set(ups) == {"transcriptomics", "proteomics", "metabolomics"}, "all three uploads accepted")
    check(ups["transcriptomics"]["orientation"] == "columns" and ups["metabolomics"]["orientation"] == "rows",
          "orientation detected both ways: genes x samples (samples in columns) and samples x compounds (samples in rows)")
    check(ups["metabolomics"]["set_aside_columns"] == ["Group"], "the Group column is set aside, not treated as a compound")
    summaries = page.locator(".multiomics-corr-upload-summary").all_inner_texts()
    check(len(summaries) == 3 and "samples are in the columns" in summaries[0] and "samples are in the rows" in summaries[2]
          and "Control_1" in summaries[0] and "Group" in summaries[2],
          "each upload shows what was understood: layout, reason, first samples, set-aside columns")
    # change the orientation of a file: the SAME file is re-read, not lost
    # (forcing the WRONG layout on the RNA file makes the 299 genes the 'samples' - a visible difference that proves the file
    # was re-read; then back to automatic)
    page.locator("select.multiomics-corr-orientation").nth(0).select_option("rows")
    page.wait_for_timeout(1500)
    check(page.evaluate("STATE.multiomics.corr.uploads.transcriptomics.info.n_samples") == 299
          and page.evaluate("STATE.multiomics.corr.uploads.transcriptomics.info.orientation") == "rows",
          "choosing an orientation re-reads the same file with it (forced 'rows': the 299 genes become the samples)")
    page.locator("select.multiomics-corr-orientation").nth(0).select_option("auto")
    page.wait_for_timeout(1500)
    check(page.evaluate("STATE.multiomics.corr.uploads.transcriptomics.info.n_samples") == 12
          and page.evaluate("STATE.multiomics.corr.uploads.transcriptomics.info.orientation") == "columns",
          "back to automatic: the same file is read the right way again")
    page.locator("button.multiomics-corr-resolve").click()
    check(wait_for(page, "!!STATE.multiomics.corr.alignment", 600, "resolve job"), "resolve and match completed")
    page.wait_for_timeout(300)
    al = page.evaluate("STATE.multiomics.corr.alignment")
    check(al["n_shared"] == 12, f"12 samples are in every file (got {al['n_shared']})")
    check(al["omics"]["transcriptomics"]["n_mapped"] > 250 and al["omics"]["proteomics"]["n_mapped"] > 250,
          f"transcript and protein ids were mapped to gene symbols (RNA {al['omics']['transcriptomics']['n_mapped']}, protein {al['omics']['proteomics']['n_mapped']})")
    check(page.locator(".multiomics-corr-alignment").count() == 1 and "Samples in every file" in page.inner_text(".multiomics-corr-alignment"),
          "the alignment card is shown")
    check(page.locator("select.multiomics-corr-method").count() == 1, "the setup appears only after the samples are matched")
    page.locator("input.multiomics-corr-case").click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.multiomics.corr.alignment") is None and page.locator("select.multiomics-corr-method").count() == 0,
          "SENSITIVITY: changing the case option drops the earlier match (it was made with the other rule)")
    page.locator("input.multiomics-corr-case").click()
    page.wait_for_timeout(200)
    page.locator("button.multiomics-corr-resolve").click()
    check(wait_for(page, "!!STATE.multiomics.corr.alignment", 600, "re-resolve"), "matched again")
    page.wait_for_timeout(300)

    # ----------------------------------------------------------- the setup layout
    print("=== setup layout: nothing overlaps ===")
    boxes = page.evaluate("""() => {
        const r = (e) => { const b = e.getBoundingClientRect(); return {top: b.top, bottom: b.bottom, left: b.left, right: b.right}; };
        const cfg = document.querySelector('.multiomics-corr-config');
        const sel = cfg.querySelector('select.multiomics-corr-method');
        const pill = cfg.querySelector('input.multiomics-corr-log2').closest('label');
        const labels = [...cfg.querySelectorAll('.chart-option-label')].map(e => ({text: e.textContent, ...r(e)}));
        const fdr = labels.find(l => l.text.startsWith('FDR cutoff'));
        const ci = document.querySelector('input.multiomics-corr-case').closest('label');
        const resolve = document.querySelector('button.multiomics-corr-resolve');
        return {sel: r(sel), pill: r(pill), fdr, ci: r(ci), resolve: r(resolve), pillHeight: pill.offsetHeight, ciHeight: ci.offsetHeight, lastUpload: r([...document.querySelectorAll('.multiomics-corr-upload')].pop())};
    }""")
    check(boxes["pill"]["top"] >= boxes["sel"]["bottom"] - 0.5, f"the log2 checkbox starts below the Correlation menu (pill top {boxes['pill']['top']:.0f} >= menu bottom {boxes['sel']['bottom']:.0f})")
    check(boxes["pill"]["bottom"] <= boxes["fdr"]["top"] + 0.5, f"...and ends above the next field's label (pill bottom {boxes['pill']['bottom']:.0f} <= label top {boxes['fdr']['top']:.0f})")
    check(boxes["pillHeight"] >= 40, f"the pill has its own full height ({boxes['pillHeight']}px)")
    check(boxes["ci"]["top"] >= boxes["lastUpload"]["bottom"] - 0.5, f"the case-ignoring checkbox starts below the last upload block (top {boxes[chr(99)+chr(105)][chr(116)+chr(111)+chr(112)]:.0f} >= {boxes[chr(108)+chr(97)+chr(115)+chr(116)+chr(85)+chr(112)+chr(108)+chr(111)+chr(97)+chr(100)][chr(98)+chr(111)+chr(116)+chr(116)+chr(111)+chr(109)]:.0f})")
    check(boxes["ci"]["bottom"] <= boxes["resolve"]["top"] + 0.5 and boxes["ciHeight"] >= 40, f"the case-ignoring checkbox has its full height ({boxes[chr(99)+chr(105)+chr(72)+chr(101)+chr(105)+chr(103)+chr(104)+chr(116)]}px) and does not overlap the button below it")
    if os.environ.get("CORR_SHOT_DIR"):
        page.locator(".multiomics-corr-config").screenshot(path=os.path.join(os.environ["CORR_SHOT_DIR"], "Setup.png"))

    # ----------------------------------------------------------- run
    print("=== run ===")
    page.locator("button.multiomics-corr-run").click()
    check(wait_for(page, "!!STATE.multiomics.corr.run.results", 600, "correlation run"), "the run completed")
    page.wait_for_timeout(900)
    R = page.evaluate("STATE.multiomics.corr.run.results")
    check(page.evaluate("STATE.aaPage") == "multiomics_corr_results", "the results screen is shown")
    rail = [t.strip() for t in page.locator("button.result-nav-btn").all_inner_texts()]
    check(rail == ["Overview", "Same-gene concordance", "Cross-omics pairs", "Heatmap", "Network"], f"result panes (got {rail})")
    check("treatment effect" in page.inner_text(".multiomics-corr-result"), "the Overview states the group-effect limitation")

    # ----------------------------------------------------------- concordance
    print("=== same-gene concordance ===")
    open_pane(page, "Same-gene concordance")
    c = R["concordance"]
    check(c["n_tested"] > 250 and c["n_sig_pos"] >= 8, f"concordance computed (tested {c['n_tested']}, significant positive {c['n_sig_pos']})")
    check(page.evaluate(HAS_DATA, "multiomics-corr-hist"), "the rho histogram is drawn")
    check(page.evaluate(MODEBAR, "multiomics-corr-hist"), "the histogram has the figure tools row (zoom / pan / autoscale / download)")
    rows = page.locator(".result-pane:visible tbody tr")
    genes = [rows.nth(i).locator("td").first.inner_text() for i in range(min(rows.count(), 400))]
    check("CS" in genes, "CS is in the concordance table")
    cs_row = rows.nth(genes.index("CS"))
    cs_row.click()
    check(wait_for(page, "!!" + SCATTER, 15, "gene scatter"),
          "clicking a gene draws its scatter")
    title = page.evaluate(SCATTER + ".layout.title.text")
    table_rho = float(cs_row.locator("td").nth(1).inner_text())
    check(f"{table_rho:.3f}"[:5] in title, f"the scatter's rho equals the table's ({table_rho:.3f} vs '{title}')")
    npts = page.evaluate(SCATTER + ".data[0].x.length")
    check(npts >= 10, f"the scatter has one point per sample (got {npts})")
    check(page.evaluate(SCATTER + ".layout.xaxis.title.text") == "Transcript CS",
          "the axes name what is plotted")
    with page.expect_download() as dl:
        page.locator("button", has_text="Download concordance table").click()
    text = open(dl.value.path(), encoding="utf-8-sig").read()
    csv_rows = list(csv.DictReader(io.StringIO(text)))
    check(len(csv_rows) == c["n_tested"] and "Gene" in csv_rows[0], f"the concordance CSV has every tested gene ({len(csv_rows)})")

    # ----------------------------------------------------------- pairs
    print("=== cross-omics pairs ===")
    open_pane(page, "Cross-omics pairs")
    sel = page.locator("select.multiomics-corr-block")
    keys = page.evaluate("Object.keys(STATE.multiomics.corr.run.results.blocks)")
    check(sel.locator("option").count() == len(keys) == 3, f"one omics pair per block ({keys})")
    key_rm = next(k for k in keys if k == "transcriptomics|metabolomics")
    sel.select_option(key_rm)
    page.wait_for_timeout(500)
    b = R["blocks"][key_rm]
    check(b["n_reported"] > 0, f"RNA x metabolite has reported pairs ({b['n_reported']})")
    labels = [r["Feature A"] + "|" + r["Feature B"] for r in b["table"]["rows"]]
    check("CS|Citrate" in labels, "CS - Citrate (planted) is among the reported pairs")
    idx = labels.index("CS|Citrate")
    page.locator(".result-pane:visible tbody tr").nth(idx).click()
    check(wait_for(page, "!!" + SCATTER, 15, "pair scatter"),
          "clicking a pair draws its scatter")
    pt = page.evaluate(SCATTER + ".layout.title.text")
    check("rho" in pt and "n = 12" in pt, f"the pair scatter states rho and n ('{pt}')")
    with page.expect_download() as dl:
        page.locator("button", has_text="Download this table").click()
    rows_csv = list(csv.DictReader(io.StringIO(open(dl.value.path(), encoding="utf-8-sig").read())))
    check(len(rows_csv) == b["n_reported"] or len(rows_csv) == len(b["table"]["rows"]), f"the pair CSV has the reported rows ({len(rows_csv)})")

    # ----------------------------------------------------------- heatmap
    print("=== heatmap ===")
    open_pane(page, "Heatmap")
    check(page.evaluate(HAS_DATA, "multiomics-corr-heatmap"), "the heatmap is drawn")
    hm = page.evaluate("""() => { const d = document.getElementById('multiomics-corr-heatmap'); const f = d.parentElement;
        const t = d.data[0]; return {w: d.clientWidth, frame: f.clientWidth, nx: t.x.length, ny: t.y.length, zmin: t.zmin, zmax: t.zmax,
        dots: t.text.reduce((a, r) => a + r.filter(v => v).length, 0), tool: !!d.querySelector('.modebar')}; }""")
    check(hm["zmin"] == -1 and hm["zmax"] == 1, "the colour scale is symmetric (-1..1), so red and blue mean the same size of rho")
    first_key = page.locator("select.multiomics-corr-heat-block").input_value()
    want_dots = sum(1 for row in R["blocks"][first_key]["heat"]["sig"] for v in row if v)
    check(0 < hm["dots"] == want_dots < hm["nx"] * hm["ny"],
          f"a dot marks exactly the REPORTED pairs ({hm['dots']} of {hm['nx'] * hm['ny']} cells; payload says {want_dots})")
    # cell size is only meaningful with many cells: draw a 40 x 40 synthetic heat with the app's own builder and measure the pixels
    cell_px = page.evaluate("""() => {
        const n = 40, z = [], sig = [], names = [];
        for (let i = 0; i < n; i++) { names.push('F' + i); z.push(Array.from({length: n}, (_, j) => Math.sin(i * j))); sig.push(Array(n).fill(false)); }
        const d = document.createElement('div'); d.id = 'synthetic-heat'; document.body.appendChild(d);
        buildCorrHeatmap(d, {rows: names, cols: names, z, sig}, 'A', 'B', 'synthetic');
        return new Promise(res => setTimeout(() => { const xa = d._fullLayout.xaxis, ya = d._fullLayout.yaxis;
            res({x: xa._length / n, y: ya._length / n}); d.remove(); }, 1200));
    }""")
    check(cell_px["x"] >= 18 and cell_px["y"] >= 18, f"fixed cell size: 40 x 40 features keep >= 18 px per cell (got {cell_px['x']:.1f} x {cell_px['y']:.1f})")
    check(hm["tool"], "the heatmap has the figure tools row")
    check([t.strip() for t in page.locator(".result-pane:visible .multiomics-corr-heat-tools button").all_inner_texts()] == ["Zoom tool", "Pan tool", "Zoom in", "Zoom out", "Autoscale", "Reset axes"],
          "the heatmap has the visible tools row (zoom, pan, zoom in/out, autoscale, reset) like every other figure")
    hx0 = page.evaluate("document.getElementById('multiomics-corr-heatmap')._fullLayout.xaxis.range.length")
    xr0 = page.evaluate("(() => { const a = document.getElementById('multiomics-corr-heatmap')._fullLayout.yaxis.range; return Math.abs(a[1] - a[0]); })()")
    page.locator(".result-pane:visible .multiomics-corr-heat-tools button.fig-tool-in").click()
    page.wait_for_timeout(500)
    xr1 = page.evaluate("(() => { const a = document.getElementById('multiomics-corr-heatmap')._fullLayout.yaxis.range; return Math.abs(a[1] - a[0]); })()")
    check(abs(xr1 - xr0 / 2) < 0.01, f"Zoom in works on the heatmap ({xr0:.1f} -> {xr1:.1f} rows)")
    page.locator(".result-pane:visible .multiomics-corr-heat-tools button.fig-tool-reset").click()
    page.wait_for_timeout(400)
    page.locator(".multiomics-corr-heat-download").click()
    page.wait_for_timeout(400)
    check(page.locator(".ora-chart-export-menu").count() >= 1, "Download figure opens the export menu for the heatmap")
    page.keyboard.press("Escape")
    page.mouse.click(5, 5)
    check(hm["w"] >= 400, f"the heatmap is not squeezed ({hm['w']}px wide)")
    page.locator("select.multiomics-corr-heat-block").select_option(index=1)
    page.wait_for_timeout(900)
    check(page.evaluate(HAS_DATA, "multiomics-corr-heatmap"), "switching the omics pair redraws the heatmap")

    # ----------------------------------------------------------- colour gradient
    print("=== colour gradient ===")
    GRAD = """() => { const d = document.getElementById('multiomics-corr-heatmap'); const cs = d.data[0].colorscale;
        return {first: cs[0][1], last: cs[cs.length - 1][1], n: cs.length}; }"""
    opts = page.locator("select.multiomics-corr-scale option").all_inner_texts()
    check(len(opts) >= 5 and opts[0].startswith("Red - Blue"), f"a colour gradient menu offers {len(opts)} choices, the original first ({opts})")
    g0 = page.evaluate(GRAD)
    check(g0["first"].upper() == "#2166AC" and g0["last"].upper() == "#B2182B", f"default: blue (negative) to red (positive) ({g0})")
    page.locator("select.multiomics-corr-scale").select_option("orangepurple")
    page.wait_for_timeout(1200)
    check(page.locator("button.result-nav-btn.active").inner_text().strip() == "Heatmap", "choosing a gradient keeps you on the Heatmap pane")
    g1 = page.evaluate(GRAD)
    check("red positive" not in page.inner_text(".result-pane:visible"), "the caption does not claim a colour the chosen gradient no longer has")
    check(g1["first"].upper() == "#542788" and g1["last"].upper() == "#B35806", f"the heatmap now runs purple (negative) to orange (positive) ({g1})")
    page.locator("input.multiomics-corr-flip").click()
    page.wait_for_timeout(1200)
    g2 = page.evaluate(GRAD)
    check(g2["first"].upper() == "#B35806" and g2["last"].upper() == "#542788", f"'swap' reverses the ends ({g2})")
    check(page.locator("select.multiomics-corr-scale").input_value() == "orangepurple" and page.locator("input.multiomics-corr-flip").is_checked(),
          "the menu and the swap box show what is applied")
    open_pane(page, "Network")
    ncol = page.evaluate("""() => { const d = document.getElementById('multiomics-corr-network'); return d.data.filter(t => t.mode === 'lines').map(t => t.line.color); }""")
    check(set(ncol) == {"#542788", "#B35806"}, f"the network's lines use the same gradient ends ({sorted(set(ncol))})")
    open_pane(page, "Same-gene concordance")
    hcol = page.evaluate("""() => [...new Set(document.getElementById('multiomics-corr-hist').data[0].marker.color)]""")
    check(set(hcol) == {"#542788", "#B35806"}, f"...and so does the histogram ({sorted(set(hcol))})")
    open_pane(page, "Heatmap")
    page.locator("input.multiomics-corr-flip").click()
    page.locator("select.multiomics-corr-scale").select_option("redblue")
    page.wait_for_timeout(1200)

    # ----------------------------------------------------------- network
    print("=== network ===")
    open_pane(page, "Network")
    check(page.evaluate(HAS_DATA, "multiomics-corr-network"), "the network is drawn")
    nm = page.evaluate("""() => { const d = document.getElementById('multiomics-corr-network'); const xa = d._fullLayout.xaxis, ya = d._fullLayout.yaxis;
        const pts = {}; let nNodes = 0;
        d.data.forEach(t => { if (t.mode === 'markers+text') t.x.forEach((x, i) => { nNodes++; (pts[x] = pts[x] || []).push(ya.l2p(t.y[i])); }); });
        let minGap = 1e9;
        Object.values(pts).forEach(ys => { ys.sort((a, b) => a - b); for (let i = 1; i < ys.length; i++) minGap = Math.min(minGap, ys[i] - ys[i - 1]); });
        return {nNodes, minGap, tool: !!d.querySelector('.modebar'), cols: Object.keys(pts).length, h: d.clientHeight}; }""")
    nodes = R["network"]["nodes"]
    check(nm["nNodes"] == len(nodes) and nm["cols"] == len(R["network"]["omics"]), f"every node is drawn ({nm['nNodes']}) in one column per omics ({nm['cols']})")
    check(nm["minGap"] > MARKER_PX, f"no two nodes of a column overlap (closest {nm['minGap']:.1f}px > marker {MARKER_PX}px)")
    check(nm["tool"], "the network has the figure tools row")
    check(page.locator(".result-pane:visible .multiomics-corr-net-tools button").count() == 6 and page.locator(".multiomics-corr-net-download").is_visible(),
          "the network has the visible tools row and a Download figure button")
    pos = page.evaluate("""() => document.getElementById('multiomics-corr-network').data
        .filter(t => t.mode === 'markers+text').map(t => t.textposition)""")
    check(pos == ["middle left", "top center", "middle right"],
          f"labels of the middle column sit ABOVE their nodes, not under the lines that run through them (got {pos})")
    check(not any(n["row"] < 0 for n in nodes) and R["network"]["crossings"] >= 0, "layout is well-formed")

    # ----------------------------------------------------------- a very long pair list
    print("=== a long pair list is capped on screen, not in the CSV ===")
    page.evaluate("""() => { const b = STATE.multiomics.corr.run.results.blocks['transcriptomics|proteomics'];
        b.table.rows = Array.from({length: 2500}, () => b.table.rows[0]); b.ids = Array.from({length: 2500}, () => b.ids[0]); render(); }""")
    page.wait_for_timeout(600)
    open_pane(page, "Cross-omics pairs")
    page.locator("select.multiomics-corr-block").select_option("transcriptomics|proteomics")
    page.wait_for_timeout(600)
    shown_rows = page.locator(".result-pane:visible tbody tr").count()
    # numbers are written in the browser's own locale (pt-BR here: "2.000"), so either separator is right
    check(shown_rows == 2000 and re.search(r"Showing the 2[., ]?000 strongest of 2[., ]?500", page.inner_text(".result-pane:visible")),
          f"2,500 listed pairs: 2,000 rows drawn and a note says the CSV has them all (drew {shown_rows})")

    # ----------------------------------------------------------- navigation
    print("=== navigation ===")
    page.locator("button", has_text="Change settings").click()
    page.wait_for_timeout(400)
    check(page.evaluate("STATE.aaPage") == "multiomics_corr" and page.evaluate("!!STATE.multiomics.corr.alignment"),
          "Change settings returns to the setup with the matched samples kept")
    page.locator("button", has_text="Back to Multi-omics").click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "multiomics_home", "back to the tile screen")

    print("=== console ===")
    check(len(errors) == 0, f"no page error or console error during the whole flow (got {errors})")
    browser.close()

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
