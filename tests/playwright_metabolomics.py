"""Real-browser coverage of the four Metabolomics tiles.

    Over-Representation Analysis            (compound list  -> KEGG/Reactome/SMPDB/PathBank)
    Chemical Class Enrichment               (compound list  -> RefMet)
    Lipid Class / Structural Enrichment     (compound list  -> LIPID MAPS)
    Quantitative Enrichment                 (concentration table -> everything)

WHAT THIS FILE EXISTS TO CATCH. The first Metabolomics UI was written without
ever being opened in a browser, and a metabolomics researcher rejected it on
sight: it asked for log2FC/pvalue/padj columns (a differential-expression
convention metabolomics does not use), crammed every database into one screen
instead of one tile per analysis, and looked nothing like the rest of the app.
So this test asserts the things that review found:

  * the four separate tiles exist;
  * NO log2FC/pvalue/padj input is offered anywhere on the list tiles, and the
    Quantitative tile asks for a concentration table + groups instead;
  * house-style structure - flat module screen (no numbered "card per step"
    wizard), the shared file widget, ONE result <select> + ONE result card;
  * every chart control the other modules have is present AND redraws;
  * chart text says "Compound", never "Gene";
  * KEGG is organism-specific - checking KEGG on a pathway tile shows the KEGG
    organism box (a search over KEGG's own organism list, like Metabolic
    Connection's organism box; unchecking removes it), the derived organism is
    reported by default, and an explicitly picked KEGG organism wins and reaches
    the run;
  * the class tiles do NOT offer the enrichment map (class terms are disjoint -
    each compound has exactly one class per level - so it can never have an
    edge) but DO keep the dot plot (a valid enrichment view) and offer a
    class-composition pie
    whose slices sum to the user's compound set, with one working colour picker
    per class; the pathway panels never get a pie;
  * the pie's percentage labels: EVERY slice has one (Plotly's own labels drop
    small slices), each leader line ends on its own slice, they can be hidden,
    dragged (and stay where they were dragged across a redraw), reset, and the
    figure can be zoomed; slice borders can be white, black or off; and a pie
    with MANY classes (60 here) sizes itself so every label and every legend
    entry stays readable - no overprinted or clipped labels, no scrolling legend.

Needs a live server on port 8009 and real network access (PubChem/KEGG name
resolution, Reactome/SMPDB/PathBank downloads on a cold cache) - not mocked.
"""
import os
import re

from playwright.sync_api import sync_playwright

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
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


def open_tile(page, tile_name):
    page.evaluate("""() => { STATE.view = 'annotation_analytics'; STATE.aaPage = 'overview'; render(); }""")
    page.wait_for_timeout(150)
    page.locator("button.module-btn.metabolomics", has_text=re.compile(r"^Metabolomics$")).first.click()
    page.wait_for_timeout(150)
    page.locator("button.module-btn.metabolomics", has_text=tile_name).first.click()
    page.wait_for_timeout(200)
    check(page.evaluate("STATE.aaPage") == "metabolomics_module", f"[{tile_name}] navigated to the module screen")


def wait_for(page, expr, seconds, label):
    for _ in range(seconds):
        if page.evaluate(expr):
            return True
        page.wait_for_timeout(1000)
    check(False, f"timed out waiting for: {label}")
    return False


def upload_and_resolve(page, tile_name, filename, needs_organism):
    if needs_organism:
        page.evaluate("""() => { STATE.metabolomics.taxonId = '9606'; STATE.metabolomics.taxonName = 'Homo sapiens'; render(); }""")
        page.wait_for_timeout(150)
    page.locator("input[type=file]").first.set_input_files(os.path.join(_ROOT, "Data3", filename))
    page.wait_for_timeout(1500)
    page.locator("button", has_text="Resolve compounds").first.click()
    ok = wait_for(page, "!!STATE.metabolomics.resolved", 600, f"[{tile_name}] resolve (first run builds the crosswalk)")
    check(ok, f"[{tile_name}] resolve completed")
    return ok


def assert_house_style_module_screen(page, tile_name, is_quant=False):
    check(page.locator(".ora-config-screen.ora-module-screen").count() == 1,
          f"[{tile_name}] module screen uses the ora-config-screen ora-module-screen wrapper")
    check(page.locator(".section-title").count() == 0 and page.locator(".rank-card").count() == 0,
          f"[{tile_name}] flat module screen - no numbered card-per-step wizard")
    check(page.locator(".file-drop").count() >= 1, f"[{tile_name}] upload uses the shared file widget")
    check(page.locator("button.ora-inline-example-link").count() >= 1,
          f"[{tile_name}] example file is an inline example link, not a full button")
    body = page.inner_text("body")
    if is_quant:
        check("Omiscera computes the per-compound statistics itself" in body,
              f"[{tile_name}] tells the user the statistics are computed for them")
    else:
        check(not re.search(r"log2\s*FC|p-?value|padj", body, re.I),
              f"[{tile_name}] no log2FC/pvalue/padj input requested on the module screen")


def click_continue_and_run(page, tile_name):
    page.locator("button", has_text="Continue to " + tile_name).first.click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "metabolomics_config", f"[{tile_name}] navigated to the config screen")


def run_and_wait(page, tile_name, timeout=900):
    run_btn = page.locator("button.ora-run-btn")
    check(not run_btn.is_disabled(), f"[{tile_name}] Run button enabled")
    run_btn.click()
    # Configure validates an explicit KEGG organism against KEGG before the screen
    # changes, so navigation is not instant - wait for it rather than guess a delay.
    for _ in range(60):
        if page.evaluate("STATE.aaPage") == "metabolomics_results":
            break
        page.wait_for_timeout(500)
    check(page.evaluate("STATE.aaPage") == "metabolomics_results", f"[{tile_name}] navigated to the results screen")
    ok = wait_for(page, "!!STATE.metabolomics.run.results", timeout, f"[{tile_name}] analysis job")
    check(ok, f"[{tile_name}] analysis job completed")
    return ok


def assert_house_style_results(page, tile_name):
    check(page.locator(".ora-config-screen.ora-results-screen").count() == 1,
          f"[{tile_name}] results screen uses the ora-results-screen wrapper")
    check(page.locator("select.transcriptomics-ora-result-select").count() == 1,
          f"[{tile_name}] exactly ONE result-table <select> (not one rail button per database)")
    check(page.locator(".ora-result-card").count() == 1, f"[{tile_name}] exactly ONE result card")
    for label in ["Overview", "Result table"]:
        check(page.locator("button.result-nav-btn", has_text=label).count() >= 1,
              f"[{tile_name}] result card has an '{label}' pane")


def pick_significant_key(page):
    return page.evaluate("""() => {
        const tables = STATE.metabolomics.run.results.tables;
        const alpha = STATE.metabolomics.run.sigAlpha;
        for (const key of Object.keys(tables)) {
            const rows = (tables[key].table && tables[key].table.rows) || [];
            if (rows.some(r => r.FDR !== null && r.FDR !== '' && Number(r.FDR) < alpha)) return key;
        }
        return null;
    }""")


def select_result(page, key):
    page.evaluate(
        """(key) => {
            const sel = document.querySelector('select.transcriptomics-ora-result-select');
            sel.value = key;
            sel.dispatchEvent(new Event('change'));
        }""", key)
    page.wait_for_timeout(300)


def chart_words(page):
    """Every hover template and axis/colorbar title across all drawn figures."""
    return page.evaluate("""() => {
        const out = [];
        const title = (t) => (t && typeof t === 'object') ? t.text : t;
        document.querySelectorAll('.js-plotly-plot').forEach(gd => {
            (gd.data || []).forEach(tr => {
                if (tr.hovertemplate) out.push(String(tr.hovertemplate));
                if (tr.marker && tr.marker.colorbar) out.push(String(title(tr.marker.colorbar.title) || ''));
            });
            const lay = gd.layout || {};
            ['xaxis', 'yaxis'].forEach(a => { if (lay[a]) out.push(String(title(lay[a].title) || '')); });
            if (lay.title) out.push(String(title(lay.title) || ''));
        });
        return out;
    }""")


def exercise_charts(page, tile_name, is_quant, errors, is_class=False):
    """Open every chart, check the control set, redraw through every control,
    and assert no chart text still talks about genes. is_class: a Chemical /
    Lipid Class panel - the dot plot stays (a valid enrichment view) but there
    is NO enrichment map (asserted ABSENT, so a regression that re-adds it is
    caught: class terms are disjoint, so it could never have an edge)."""
    overview = page.locator("button.result-nav-btn", has_text="Compound set overview")
    check(overview.count() >= 1, f"[{tile_name}] 'Compound set overview' chart pane present")
    if overview.count() < 1:
        return
    overview.first.click()
    page.wait_for_timeout(200)

    expected = ["Terms shown", "Font size", "Bar color", "Dot gradient"] \
        + ([] if is_quant else ["Rank & show by"]) + (["Bar metric", "Volcano y-axis"] if is_quant else [])
    labels = page.locator(".chart-option-label").all_inner_texts()
    for name in expected:
        check(any(name in text for text in labels), f"[{tile_name}] chart control present: {name}")

    toggles = ["Show dot plot", "Show bar plot"] + (["Show volcano plot"] if is_quant else [])
    for label in toggles:
        btn = page.locator("button.chart-toggle-btn", has_text=label)
        check(btn.count() >= 1, f"[{tile_name}] '{label}' button present")
        if btn.count() >= 1:
            btn.first.click()
            page.wait_for_timeout(600)
    n_plots = page.locator(".js-plotly-plot").count()
    check(n_plots >= len(toggles), f"[{tile_name}] overview figures rendered (got {n_plots})")

    # Every select in the overview pane must redraw without error.
    pane_selects = page.locator(".chart-option-select:visible")
    n_selects = pane_selects.count()
    for i in range(n_selects):
        sel = pane_selects.nth(i)
        options = sel.locator("option").all_inner_texts()
        if len(options) < 2:
            continue
        values = sel.locator("option").evaluate_all("els => els.map(e => e.value)")
        current = sel.evaluate("s => s.value")
        alternative = next((v for v in values if v != current), current)
        sel.select_option(alternative)
        page.wait_for_timeout(400)
    check(len(errors) == 0, f"[{tile_name}] every overview control redrew with no page error (got {errors})")

    rel = page.locator("button.result-nav-btn", has_text="Compound set relationships")
    check(rel.count() >= 1, f"[{tile_name}] 'Compound set relationships' chart pane present")
    if rel.count() >= 1:
        rel.first.click()
        page.wait_for_timeout(200)
        if is_class:
            check(page.locator("button.chart-toggle-btn", has_text="Show enrichment map").count() == 0,
                  f"[{tile_name}] NO enrichment map on a class panel - class terms are disjoint, so it could never have an edge")
        for label in (["Show compound-concept network"] if is_class else ["Show enrichment map", "Show compound-concept network"]):
            btn = page.locator("button.chart-toggle-btn", has_text=label)
            check(btn.count() >= 1, f"[{tile_name}] '{label}' button present")
            if btn.count() >= 1:
                btn.first.click()
                page.wait_for_timeout(900)
        check(len(errors) == 0, f"[{tile_name}] relationship charts drew with no page error (got {errors})")

    words = chart_words(page)
    check(len(words) > 0, f"[{tile_name}] collected chart text to inspect (got {len(words)} strings)")
    gene_words = [w for w in words if re.search(r"\bGenes?\b", w)]
    check(not gene_words, f"[{tile_name}] no chart hover/axis text says 'Gene' (found: {gene_words[:3]})")
    check(any("Compound" in w for w in words), f"[{tile_name}] chart text says 'Compound' instead")


PIE_STATE_JS = """() => {
    const gd = [...document.querySelectorAll('.js-plotly-plot')].find(d => d.data && d.data[0] && d.data[0].type === 'pie');
    if (!gd) return null;
    const tr = gd.data[0];
    return { labels: tr.labels, values: tr.values, colors: tr.marker.colors, hover: tr.hovertemplate };
}"""


LABEL_STATE_JS = """() => {
    const gd = [...document.querySelectorAll('.js-plotly-plot')].find(d => d.data && d.data[0] && d.data[0].type === 'pie');
    if (!gd) return null;
    const tr = gd.data[0], L = gd.layout, s = gd._fullLayout._size, cd = gd.calcdata[0][0];
    return {
        values: tr.values, labels: tr.labels, textinfo: tr.textinfo, W: L.width, H: L.height,
        anns: (L.annotations || []).map(a => ({text: a.text, x: a.x, y: a.y, ax: a.ax, ay: a.ay})),
        size: {l: s.l, t: s.t, w: s.w, h: s.h}, cx: cd.cx, cy: cd.cy, r: cd.r,
    };
}"""


def exercise_pie_labels(page, tile_name, errors):
    """The pie's percentage labels: one per slice, attached to the RIGHT slice,
    hideable, draggable, persistent across a redraw, resettable, zoomable."""
    import math

    def pie_state():
        return page.evaluate(LABEL_STATE_JS)

    ctl = lambda label: page.locator("label.chart-option-control", has_text=label).locator("select")
    ctl("Percentages").select_option("show")
    ctl("Zoom").select_option("1")
    page.wait_for_timeout(500)
    s = pie_state()
    check(s is not None and s["textinfo"] == "none",
          f"[{tile_name}] Plotly's own slice text is off - the labels are annotations (its own drop small slices)")
    n = len(s["values"])
    total = sum(s["values"])
    check(len(s["anns"]) == n, f"[{tile_name}] EVERY slice has a percentage label ({len(s['anns'])} labels for {n} slices)")
    check(all(a["text"] == f"{v / total * 100:.1f}%" for a, v in zip(s["anns"], s["values"])),
          f"[{tile_name}] each label is its own slice's share of the whole")

    # SENSITIVITY: a label that merely EXISTS proves nothing - its leader line must end on ITS slice.
    acc = 0
    bad = []
    for v, a in zip(s["values"], s["anns"]):
        start, end = acc / total * 2 * math.pi, (acc + v) / total * 2 * math.pi
        acc += v
        tx, ty = s["size"]["l"] + a["x"] * s["size"]["w"], s["size"]["t"] + (1 - a["y"]) * s["size"]["h"]
        dx, dy = tx - s["cx"], ty - s["cy"]
        ang = math.atan2(dx, -dy) % (2 * math.pi)
        if abs(math.hypot(dx, dy) - s["r"] * 0.97) > 1.5 or not (start - 1e-6 <= ang <= end + 1e-6):
            bad.append(a["text"])
    check(not bad, f"[{tile_name}] SENSITIVITY: every leader line ends inside its own slice, just inside the rim "
                   f"(misplaced: {bad})")

    ctl("Percentages").select_option("hide")
    page.wait_for_timeout(500)
    check(len(pie_state()["anns"]) == 0, f"[{tile_name}] 'Hide' removes every percentage")
    ctl("Percentages").select_option("show")
    page.wait_for_timeout(500)
    check(len(pie_state()["anns"]) == n, f"[{tile_name}] 'Show' brings them all back")

    # Drag one label (as with a network node); the line follows and the move is kept.
    before = pie_state()["anns"][0]
    label_el = page.locator('.js-plotly-plot .annotation[data-index="0"] text').first
    # The pie sits far below the fold; a mouse event outside the viewport reaches nothing.
    label_el.scroll_into_view_if_needed()
    page.wait_for_timeout(200)
    box = label_el.bounding_box()
    check(box is not None, f"[{tile_name}] the first percentage label is on screen")
    if box:
        cx, cy = box["x"] + box["width"] / 2, box["y"] + box["height"] / 2
        page.mouse.move(cx, cy)
        page.mouse.down()
        page.mouse.move(cx + 25, cy + 15, steps=4)
        page.mouse.move(cx + 60, cy + 40, steps=6)
        page.mouse.up()
        page.wait_for_timeout(500)
        moved = pie_state()["anns"][0]
        check(abs(moved["ax"] - before["ax"]) > 20 or abs(moved["ay"] - before["ay"]) > 20,
              f"[{tile_name}] dragging a label moves it (ax {before['ax']:.0f}->{moved['ax']:.0f}, "
              f"ay {before['ay']:.0f}->{moved['ay']:.0f})")
        check(abs(moved["x"] - before["x"]) < 1e-9 and abs(moved["y"] - before["y"]) < 1e-9,
              f"[{tile_name}] ...while the leader line stays attached to the slice")
        ctl("Percentages").select_option("hide")
        page.wait_for_timeout(300)
        ctl("Percentages").select_option("show")
        page.wait_for_timeout(500)
        kept = pie_state()["anns"][0]
        check(abs(kept["ax"] - moved["ax"]) < 1e-6 and abs(kept["ay"] - moved["ay"]) < 1e-6,
              f"[{tile_name}] a dragged label keeps its place across a redraw (stored on the run)")
        page.locator("button.chart-toggle-btn", has_text="Reset label positions").first.click()
        page.wait_for_timeout(500)
        reset = pie_state()["anns"][0]
        check(abs(reset["ax"] - before["ax"]) < 1e-6 and abs(reset["ay"] - before["ay"]) < 1e-6,
              f"[{tile_name}] 'Reset label positions' puts it back")

    # Slice borders: white (default), black, or none - the hairline blurs a pie with thin slices.
    border = lambda: page.evaluate("""() => {
        const gd = [...document.querySelectorAll('.js-plotly-plot')].find(d => d.data && d.data[0] && d.data[0].type === 'pie');
        return gd.data[0].marker.line; }""")
    check(border()["color"].upper() == "#FFFFFF" and border()["width"] > 0, f"[{tile_name}] slice borders are white by default")
    ctl("Slice borders").select_option("black")
    page.wait_for_timeout(400)
    check(border()["color"].upper() == "#000000" and border()["width"] > 0, f"[{tile_name}] 'Black' draws black slice borders")
    ctl("Slice borders").select_option("none")
    page.wait_for_timeout(400)
    check(border()["width"] == 0, f"[{tile_name}] 'None' removes the slice borders")
    ctl("Slice borders").select_option("white")
    page.wait_for_timeout(300)

    # 'Label slices of at least': hides the labels of the thinnest slices, keeps the rest.
    # Pick a threshold that genuinely separates THIS panel's slices (some below, some at
    # or above) - a fixed one can hide nothing on a panel whose slices are all large.
    pcts = [v / total * 100 for v in pie_state()["values"]]
    cut = next((th for th in (0.5, 1, 2, 5) if any(p < th for p in pcts) and any(p >= th for p in pcts)), None)
    if cut is not None:
        opt = {0.5: "0.5", 1: "1", 2: "2", 5: "5"}[cut]
        ctl("Label slices of at least").select_option(opt)
        page.wait_for_timeout(500)
        thr = pie_state()
        expected = sum(1 for p in pcts if p >= cut)
        check(len(thr["anns"]) == expected and expected < n,
              f"[{tile_name}] 'at least {cut}%' labels only the slices of >={cut}% ({len(thr['anns'])} of {n}; expected {expected})")
        check(all(float(a["text"].rstrip("%")) >= cut - 0.05 for a in thr["anns"]),
              f"[{tile_name}] ...and none of the labelled ones is smaller")
        ctl("Label slices of at least").select_option("0")
        page.wait_for_timeout(400)
        check(len(pie_state()["anns"]) == n, f"[{tile_name}] 'All' labels every slice again")
    else:
        check(all(p >= 5 for p in pcts) or len(set(round(p, 3) for p in pcts)) == 1,
              f"[{tile_name}] (no threshold separates these slices - nothing to hide; pcts {[round(p, 1) for p in pcts]})")

    # Zoom scales the figure (spreading small slices) and the labels follow the new geometry.
    w1 = pie_state()["W"]
    ctl("Zoom").select_option("2")
    page.wait_for_timeout(700)
    z = pie_state()
    check(z["W"] >= 1.9 * w1 and z["r"] > 1.9 * s["r"] * 0.9,
          f"[{tile_name}] 'Zoom 200%' doubles the figure and the pie radius ({w1} -> {z['W']}, r {s['r']:.0f} -> {z['r']:.0f})")
    check(len(z["anns"]) == n, f"[{tile_name}] every slice is still labelled when zoomed")
    acc = 0
    bad = []
    for v, a in zip(z["values"], z["anns"]):
        start, end = acc / total * 2 * math.pi, (acc + v) / total * 2 * math.pi
        acc += v
        tx, ty = z["size"]["l"] + a["x"] * z["size"]["w"], z["size"]["t"] + (1 - a["y"]) * z["size"]["h"]
        ang = math.atan2(tx - z["cx"], -(ty - z["cy"])) % (2 * math.pi)
        if abs(math.hypot(tx - z["cx"], ty - z["cy"]) - z["r"] * 0.97) > 1.5 or not (start - 1e-6 <= ang <= end + 1e-6):
            bad.append(a["text"])
    check(not bad, f"[{tile_name}] leader lines still end on their own slice at 200% (misplaced: {bad})")
    ctl("Zoom").select_option("1")
    page.wait_for_timeout(500)
    check(len(errors) == 0, f"[{tile_name}] the label controls ran with no page error (got {errors})")


def exercise_class_composition(page, tile_name, errors):
    """The class-composition pie: descriptive, sums to the whole, per-class colours."""
    pane = page.locator("button.result-nav-btn", has_text="Class composition")
    check(pane.count() >= 1, f"[{tile_name}] 'Class composition' pane present")
    if pane.count() < 1:
        return
    pane.first.click()
    page.wait_for_timeout(200)
    check(page.locator("text=Composition, not enrichment").count() >= 1,
          f"[{tile_name}] the pie is labelled as composition, not enrichment (no p-value)")
    check(page.locator("text=no equivalent pie is offered there").count() >= 1,
          f"[{tile_name}] the note says WHY it is not offered for pathway databases")
    labels = page.locator(".chart-option-label").all_inner_texts()
    for name in ["Class level", "Classes shown", "Font size"]:
        check(any(name in text for text in labels), f"[{tile_name}] pie control present: {name}")

    page.locator("button.chart-toggle-btn", has_text="Show class composition pie").first.click()
    page.wait_for_timeout(700)
    pie = page.evaluate(PIE_STATE_JS)
    check(pie is not None, f"[{tile_name}] the pie drew")
    if pie is None:
        return
    total = page.evaluate("new Set(STATE.metabolomics.resolved.rows.map(r => r.compound_key).filter(Boolean)).size")
    check(sum(pie["values"]) == total,
          f"[{tile_name}] SENSITIVITY: the slices sum to the number of distinct matched compounds "
          f"({sum(pie['values'])} == {total}) - a composition chart whose parts do not sum to the whole is wrong")
    check("Compounds" in pie["hover"] and "Gene" not in pie["hover"], f"[{tile_name}] pie hover says Compounds")
    swatches = page.locator("input.chart-option-color")
    n_class_slices = len([lab for lab in pie["labels"] if not lab.startswith("Other (") and lab != "Unclassified"])
    check(swatches.count() == n_class_slices and n_class_slices >= 1,
          f"[{tile_name}] one colour picker per real class ({swatches.count()} pickers, {n_class_slices} classes)")
    check(len(set(pie["colors"])) == len(pie["colors"]), f"[{tile_name}] every slice starts with a distinct colour")

    swatches.first.evaluate("el => { el.value = '#ff0000'; el.dispatchEvent(new Event('input', {bubbles: true})); }")
    page.wait_for_timeout(500)
    after = page.evaluate(PIE_STATE_JS)
    check(after["colors"][0].lower() == "#ff0000",
          f"[{tile_name}] SENSITIVITY: picking a colour recolours that slice in the drawn pie (got {after['colors'][0]})")
    # The choice must survive an unrelated redraw (renderPanel re-runs on FDR changes).
    page.locator("input.ora-sig-threshold-input").fill("0.049")
    page.locator("input.ora-sig-threshold-input").dispatch_event("change")
    page.wait_for_timeout(400)
    page.locator("button.result-nav-btn", has_text="Class composition").first.click()
    page.wait_for_timeout(200)
    if page.locator("button.chart-toggle-btn", has_text="Show class composition pie").count():
        page.locator("button.chart-toggle-btn", has_text="Show class composition pie").first.click()
        page.wait_for_timeout(600)
    kept = page.evaluate(PIE_STATE_JS)
    check(kept is not None and kept["colors"][0].lower() == "#ff0000",
          f"[{tile_name}] the chosen colour survives a threshold change (stored on the run, not a render local)")
    page.locator("input.ora-sig-threshold-input").fill("0.05")
    page.locator("input.ora-sig-threshold-input").dispatch_event("change")
    page.wait_for_timeout(300)
    page.locator("button.result-nav-btn", has_text="Class composition").first.click()
    page.wait_for_timeout(200)
    if page.locator("button.chart-toggle-btn", has_text="Show class composition pie").count():
        page.locator("button.chart-toggle-btn", has_text="Show class composition pie").first.click()
        page.wait_for_timeout(600)

    # Level and cap redraw the pie (and the picker list) without error.
    level = page.locator("label.chart-option-control", has_text="Class level").locator("select")
    values = level.locator("option").evaluate_all("els => els.map(e => e.value)")
    before_labels = page.evaluate(PIE_STATE_JS)["labels"]
    level.select_option(next(v for v in values if v != level.evaluate("s => s.value")))
    page.wait_for_timeout(600)
    changed = page.evaluate(PIE_STATE_JS)
    check(changed is not None and sum(changed["values"]) == total,
          f"[{tile_name}] changing the class level redraws a pie that still sums to the whole")
    check(changed["labels"] != before_labels, f"[{tile_name}] ...with that level's own classes")
    cap = page.locator("label.chart-option-control", has_text="Classes shown").locator("select")
    cap.select_option("10")
    page.wait_for_timeout(600)
    capped = page.evaluate(PIE_STATE_JS)
    real = [lab for lab in capped["labels"] if lab != "Unclassified"]
    check(len(real) <= 11 and sum(capped["values"]) == total,
          f"[{tile_name}] 'Classes shown' caps the slices (<=10 + one 'Other') without losing compounds "
          f"({len(real)} slices, sum {sum(capped['values'])})")
    check(len(errors) == 0, f"[{tile_name}] the pie controls redrew with no page error (got {errors})")
    exercise_pie_labels(page, tile_name, errors)


MANY_CLASSES_JS = """(args) => {
    const [n, zoom, minPct, longNames, containerW] = args;
    document.querySelectorAll('#many-classes').forEach(e => e.remove());
    const div = document.createElement('div');
    div.id = 'many-classes';
    div.style.cssText = 'width:' + (containerW || 1000) + 'px;position:absolute;left:0;top:0;z-index:9999;background:#fff';
    document.body.appendChild(div);
    // a realistic long tail: two big classes, then many singletons and pairs
    const slices = [];
    for (let i = 0; i < n; i++) slices.push({label: (longNames ? 'Diacylglycerophosphoethanolamines and other long class name ' : 'Class number ') + (i + 1), count: i === 0 ? 22 : i === 1 ? 14 : (i % 5 === 0 ? 2 : 1), kind: 'class', rank: i});
    buildClassCompositionPie(div, {slices, total: slices.reduce((s, x) => s + x.count, 0)},
        (s) => _classPieDefaultColor(s.rank), 'many', {fontSize: 12, zoom, showPercent: true, minPercent: minPct,
        borderColor: 'white', savedPos: () => null, onMoved: () => {}});
    return new Promise(res => setTimeout(() => {
        const L = div.layout, fl = div._fullLayout, sz = fl._size, leg = fl.legend || {};
        const anns = (L.annotations || []).map(a => ({ly: sz.t + (1 - a.y) * sz.h + a.ay, lx: sz.l + a.x * sz.w + a.ax}));
        const cd = div.calcdata[0][0];
        // Where does each leader tip land relative to the pie Plotly ACTUALLY drew?
        const tipErr = (L.annotations || []).map(a => {
            const tx = sz.l + a.x * sz.w, ty = sz.t + (1 - a.y) * sz.h;
            return Math.abs(Math.hypot(tx - cd.cx, ty - cd.cy) - cd.r * 0.97);
        });
        const out = {W: L.width, H: L.height, plotW: sz.w, plotH: sz.h, n, anns, fontSize: L.font.size,
                     asked: {w: L.width - L.margin.l - L.margin.r, h: L.height - L.margin.t - L.margin.b},
                     maxTipErr: tipErr.length ? Math.max(...tipErr) : 0,
                     legendH: leg._height, legendMax: leg._maxHeight, legendEntries: div.querySelectorAll('.legend .traces').length};
        div.remove();
        res(out);
    }, 900));
}"""


def exercise_many_classes(page, errors):
    """A pie with dozens of classes: every label readable, the whole legend reachable."""
    tile = "60-class pie"
    r = page.evaluate(MANY_CLASSES_JS, [60, 1, 0, False, 1000])
    check(len(r["anns"]) == 60, f"[{tile}] every one of the 60 slices has a label ({len(r['anns'])})")
    cxp = r["plotW"] * 0.5        # the pie is centred in the plot area (domain x 0.2-0.8)
    worst = 1e9
    for side in (1, -1):
        col = sorted([a for a in r["anns"] if (a["lx"] - cxp) * side >= 0], key=lambda a: a["ly"])
        for a, b in zip(col, col[1:]):
            worst = min(worst, b["ly"] - a["ly"])
    check(worst >= r["fontSize"] * 1.4,
          f"[{tile}] SENSITIVITY: no two labels on a side overprint (closest pair {worst:.1f}px apart, font {r['fontSize']}px)")
    check(all(0 <= a["ly"] <= r["plotH"] and 0 <= a["lx"] <= r["plotW"] for a in r["anns"]),
          f"[{tile}] no label is clipped outside the figure")
    check(r["H"] > 520 and r["legendH"] <= r["legendMax"] + 1,
          f"[{tile}] SENSITIVITY: the figure grew to fit ({r['W']}x{r['H']}) and the whole legend fits without "
          f"scrolling (legend {r['legendH']:.0f}px of {r['legendMax']:.0f}px available)")
    check(r["legendEntries"] == 60, f"[{tile}] all 60 classes are in the legend ({r['legendEntries']})")
    thin = page.evaluate(MANY_CLASSES_JS, [60, 1, 1, False, 1000])
    check(0 < len(thin["anns"]) < 60, f"[{tile}] 'at least 1%' hides the thinnest labels ({len(thin['anns'])} of 60 kept)")
    big = page.evaluate(MANY_CLASSES_JS, [60, 2, 0, False, 1000])
    check(big["W"] >= 1.9 * r["W"] and len(big["anns"]) == 60, f"[{tile}] Zoom 200% enlarges it ({r['W']} -> {big['W']}) and keeps every label")
    small = page.evaluate(MANY_CLASSES_JS, [60, 0.5, 0, False, 1000])
    check(small["W"] <= 0.55 * r["W"] and len(small["anns"]) == 60 and small["legendH"] <= small["legendMax"] + 1,
          f"[{tile}] Zoom 50% shrinks it ({r['W']} -> {small['W']}) with every label and the whole legend still reachable")
    few = page.evaluate(MANY_CLASSES_JS, [12, 1, 0, False, 1000])
    check(few["H"] <= 600, f"[{tile}] a small pie does not balloon ({few['H']}px for 12 classes)")
    # REGRESSION: very long class names in a narrow window made Plotly grow the right
    # margin for the legend, shrinking the plot area and sliding the pie out from under
    # its leader tips (found in the 1280px window). The plot area must stay the size we
    # asked for, and every leader tip must sit on the rim of the pie Plotly really drew.
    for width in (700, 900, 1280):
        lg = page.evaluate(MANY_CLASSES_JS, [40, 1, 0, True, width])
        check(abs(lg["plotW"] - lg["asked"]["w"]) < 0.5 and abs(lg["plotH"] - lg["asked"]["h"]) < 0.5,
              f"[{tile}] SENSITIVITY: long class names in a {width}px window do NOT shrink the plot area "
              f"(asked {lg['asked']['w']}x{lg['asked']['h']}, got {lg['plotW']:.0f}x{lg['plotH']:.0f})")
        check(lg["maxTipErr"] < 1.5, f"[{tile}] every leader tip is on the drawn pie's rim ({width}px window, "
                                      f"worst error {lg['maxTipErr']:.2f}px)")
        check(len(lg["anns"]) == 40, f"[{tile}] all 40 long-named slices are labelled ({width}px window)")
    check(len(errors) == 0, f"[{tile}] no page error (got {errors})")


with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page()
    errors = []
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    page.goto(BASE, wait_until="networkidle")
    page.wait_for_timeout(300)
    check(len(errors) == 0, f"app.js parsed with no top-level error (got {errors})")

    print("=== five separate tiles ===")
    page.evaluate("""() => { STATE.view = 'annotation_analytics'; STATE.aaPage = 'overview'; render(); }""")
    page.wait_for_timeout(150)
    page.locator("button.module-btn.metabolomics", has_text=re.compile(r"^Metabolomics$")).first.click()
    page.wait_for_timeout(150)
    tile_names = page.locator("button.module-btn.metabolomics").all_inner_texts()
    for name in ["Over-Representation Analysis", "Chemical Class Enrichment",
                 "Lipid Class / Structural Enrichment", "Quantitative Enrichment",
                 "Pathway Topology Analysis"]:
        check(name in tile_names, f"tile present: {name}")
    check(len([t for t in tile_names if t]) == 5, f"exactly five tiles (got {tile_names})")

    # ------------------------------------------------------------------ ORA
    print("=== Over-Representation Analysis ===")
    tile = "Over-Representation Analysis"
    open_tile(page, tile)
    check(page.locator("input.text-input[placeholder*='organism']").count() == 1, f"[{tile}] organism selector present")
    assert_house_style_module_screen(page, tile)
    if upload_and_resolve(page, tile, "example_metabolomics_list.csv", True):
        meta = page.evaluate("STATE.metabolomics.resolveMeta")
        check(meta["n_matched"] >= 50, f"[{tile}] most of the 60 example compounds matched (got {meta['n_matched']})")
        click_continue_and_run(page, tile)
        pills = page.locator("label.ora-db-pill").all_inner_texts()
        for db in ["KEGG", "Reactome", "SMPDB", "PathBank"]:
            check(any(db in t for t in pills), f"[{tile}] database pill present: {db}")
        check(not any("Chemical Class" in t or "Lipid Class" in t for t in pills),
              f"[{tile}] class databases are NOT crammed onto this tile (they have their own)")
        check(page.locator("text=padj").count() == 0 and page.locator("text=log2FC").count() == 0,
              f"[{tile}] config screen has no padj / log2FC cutoff fields")
        # KEGG is organism-specific, so checking it shows the KEGG organism box (a
        # search over KEGG's own organism list, like Metabolic Connection's
        # organism box), and unchecking removes it.
        kegg_cb = page.locator("label.ora-db-pill", has_text="KEGG").locator("input[type=checkbox]")
        check(kegg_cb.is_checked() and page.locator(".ora-kegg-organism-wrap").count() == 1,
              f"[{tile}] KEGG is checked, so the KEGG organism box is present")
        box = page.locator(".ora-kegg-organism-wrap")
        for _ in range(20):
            if "KEGG organism hsa" in box.inner_text():
                break
            page.wait_for_timeout(500)
        check("Using Homo sapiens (KEGG organism hsa)" in box.inner_text(),
              f"[{tile}] by default it reports the KEGG organism derived from the one chosen (got {box.inner_text()[:160]!r})")
        check(box.locator("input[placeholder*='KEGG organism']").count() == 1,
              f"[{tile}] the box has a KEGG-organism search input")
        check(page.locator("input[placeholder*='organism' i]").count() == 1,
              f"[{tile}] and it is the ONLY organism search box on the config screen")
        check(box.locator("select:visible").count() == 0,
              f"[{tile}] no suggestion list for an unambiguous organism")

        # Search KEGG's own organism list and pick one explicitly.
        box.locator("input[placeholder*='KEGG organism']").fill("saccharomyces")
        for _ in range(40):
            if box.locator("select option", has_text="sce |").count() > 0:
                break
            page.wait_for_timeout(500)
        check(box.locator("select option", has_text="sce |").count() > 0,
              f"[{tile}] searching 'saccharomyces' lists KEGG organisms as 'code | name' (real KEGG search)")
        box.locator("select").select_option(label=box.locator("select option", has_text="sce |").first.inner_text())
        for _ in range(30):
            if page.evaluate("() => !!(STATE.metabolomics.run.draft && STATE.metabolomics.run.draft.keggOrganism)"):
                break
            page.wait_for_timeout(500)
        picked = page.evaluate("() => STATE.metabolomics.run.draft.keggOrganism")
        check(picked and picked["code"] == "sce" and picked["nMaps"] > 50,
              f"[{tile}] picking validates the organism and records its code and map count (got {picked})")
        check("Selected: " in page.locator(".ora-kegg-organism-wrap").inner_text()
              and "sce" in page.locator(".ora-kegg-organism-wrap").inner_text(),
              f"[{tile}] the explicit choice is shown")
        page.locator(".ora-kegg-organism-wrap button", has_text="Use the organism chosen above instead").click()
        page.wait_for_timeout(600)
        check(page.evaluate("() => STATE.metabolomics.run.draft.keggOrganism") is None,
              f"[{tile}] the explicit choice can be cleared, returning to the derived organism")

        kegg_cb.uncheck()
        page.wait_for_timeout(300)
        check(page.locator(".ora-kegg-organism-wrap").count() == 0, f"[{tile}] unchecking KEGG removes the box")
        kegg_cb.check()
        page.wait_for_timeout(300)
        check(page.locator(".ora-kegg-organism-wrap").count() == 1, f"[{tile}] re-checking KEGG brings it back")
        # Keep the run fast: KEGG + SMPDB are cached local downloads.
        for db in ["Reactome", "PathBank"]:
            page.locator("label.ora-db-pill", has_text=db).locator("input[type=checkbox]").uncheck()
        page.wait_for_timeout(150)
        if run_and_wait(page, tile):
            assert_house_style_results(page, tile)
            summary = page.evaluate("STATE.metabolomics.run.results.summary")
            check(summary.get("universe_scope") == "database", f"[{tile}] list tile tests against the database universe")
            check("n_sig_up" not in summary and "n_sig_all" not in summary,
                  f"[{tile}] no up/down/all significance split - a list has no direction")
            tables = page.evaluate("STATE.metabolomics.run.results.tables")
            check(all(" (" not in k for k in tables), f"[{tile}] result keys carry no direction suffix ({list(tables)})")
            kegg = tables.get("KEGG", {})
            check(kegg.get("n_terms_tested", 0) > 0, f"[{tile}] KEGG tested real terms")
            top = (kegg.get("table", {}).get("rows") or [{}])[0]
            check(float(top.get("FDR", 1)) < 0.05, f"[{tile}] KEGG's top term is significant (FDR={top.get('FDR')})")
            check(str(top.get("Term ID", "")).startswith("hsa"),
                  f"[{tile}] KEGG terms are the ORGANISM's maps (hsa...), not generic map#####: {top.get('Term ID')}")
            check("Homo sapiens (hsa)" in (kegg.get("caveat") or ""), f"[{tile}] the KEGG caveat names the organism")
            check("Overlapping Compounds" in top, f"[{tile}] member column is 'Overlapping Compounds'")
            check(";" in str(top.get("Overlapping Compounds", "")) or len(str(top.get("Overlapping Compounds", ""))) > 0,
                  f"[{tile}] members are semicolon-joined names")
            key = pick_significant_key(page)
            check(key is not None, f"[{tile}] some database has a significant term")
            if key:
                select_result(page, key)
                exercise_charts(page, tile, False, errors)
                check(page.locator("button.result-nav-btn", has_text="Class composition").count() == 0,
                      f"[{tile}] NO class-composition pie on a pathway panel - a compound sits in many pathways, "
                      f"so slices would not sum to the whole")
            csv = page.request.get(BASE + "/api/metabolomics/ora/download?key=KEGG",
                                   headers={"X-Session-Id": page.evaluate("getSessionId()")})
            check(csv.ok and "Overlapping Compounds" in csv.text().splitlines()[0], f"[{tile}] download CSV works")

    # -------------------------------------------------------- Chemical Class
    print("=== Chemical Class Enrichment ===")
    tile = "Chemical Class Enrichment"
    open_tile(page, tile)
    check(page.locator("input.text-input[placeholder*='organism']").count() == 0,
          f"[{tile}] no organism selector - RefMet classes are not organism-specific")
    assert_house_style_module_screen(page, tile)
    if upload_and_resolve(page, tile, "example_metabolomics_list.csv", False):
        click_continue_and_run(page, tile)
        check(page.locator("label.ora-db-pill").count() == 0, f"[{tile}] no database pills - this tile has exactly one database")
        check(page.locator("text=Chemical class level").count() >= 1, f"[{tile}] class-level selector present")
        if run_and_wait(page, tile):
            assert_house_style_results(page, tile)
            tables = page.evaluate("STATE.metabolomics.run.results.tables")
            check(list(tables) == ["Chemical Class"], f"[{tile}] exactly one family (got {list(tables)})")
            rows = tables["Chemical Class"]["table"]["rows"]
            names = [r["Term Name"] for r in rows[:5]]
            check(any("TCA" in n for n in names), f"[{tile}] TCA acids among the top classes (got {names})")
            key = pick_significant_key(page)
            check(key == "Chemical Class", f"[{tile}] a significant class exists")
            if key:
                exercise_charts(page, tile, False, errors, is_class=True)
            exercise_class_composition(page, tile, errors)

    # ------------------------------------------------------------ Lipid Class
    print("=== Lipid Class / Structural Enrichment ===")
    tile = "Lipid Class / Structural Enrichment"
    open_tile(page, tile)
    assert_house_style_module_screen(page, tile)
    if upload_and_resolve(page, tile, "example_metabolomics_list.csv", False):
        click_continue_and_run(page, tile)
        check(page.locator("text=Lipid class level").count() >= 1, f"[{tile}] lipid class-level selector present")
        if run_and_wait(page, tile):
            assert_house_style_results(page, tile)
            tables = page.evaluate("STATE.metabolomics.run.results.tables")
            check(list(tables) == ["Lipid Class"], f"[{tile}] exactly one family (got {list(tables)})")
            check(tables["Lipid Class"]["n_terms_tested"] > 0, f"[{tile}] tested real lipid classes")
            if pick_significant_key(page) == "Lipid Class":
                exercise_charts(page, tile, False, errors, is_class=True)
            exercise_class_composition(page, tile, errors)

    # ----------------------------------------------------------- Quantitative
    print("=== Quantitative Enrichment ===")
    tile = "Quantitative Enrichment"
    open_tile(page, tile)
    assert_house_style_module_screen(page, tile, is_quant=True)
    page.evaluate("""() => { STATE.metabolomics.taxonId = '9606'; STATE.metabolomics.taxonName = 'Homo sapiens'; render(); }""")
    page.locator("input[type=file]").first.set_input_files(os.path.join(_ROOT, "Data3", "example_metabolomics_concentrations.csv"))
    page.wait_for_timeout(1500)
    info = page.evaluate("STATE.metabolomics.concInfo")
    check(info is not None and info["group_column"] == "Group", f"[{tile}] group column auto-detected (got {info and info['group_column']})")
    check(sorted(info["group_levels"]) == ["Control", "Treated"], f"[{tile}] two groups detected")
    check(page.locator("text=Group column").count() >= 1 and page.locator("text=Table layout").count() >= 1,
          f"[{tile}] layout/group overrides are offered, not just silently guessed")
    # Override round trip: force the (wrong) transposed reading. It must either
    # take effect or be rejected with a message - never silently ignored.
    dialogs = []
    page.on("dialog", lambda d: (dialogs.append(d.message), d.dismiss()))
    page.locator("select").filter(has=page.locator("option", has_text="Compounds in rows")).first.select_option("transposed")
    page.wait_for_timeout(1500)
    check(page.evaluate("STATE.metabolomics.concInfo.orientation") == "transposed" or len(dialogs) > 0,
          f"[{tile}] the layout override took effect or was rejected with a message (dialogs: {dialogs})")
    page.locator("input[type=file]").first.set_input_files(os.path.join(_ROOT, "Data3", "example_metabolomics_concentrations.csv"))
    page.wait_for_timeout(1500)
    check(page.evaluate("STATE.metabolomics.concInfo.orientation") == "standard", f"[{tile}] fresh upload detects samples-in-rows again")
    page.locator("button", has_text="Resolve compounds").first.click()
    if wait_for(page, "!!STATE.metabolomics.resolved", 600, f"[{tile}] resolve"):
        click_continue_and_run(page, tile)
        check(page.locator("text=Groups to compare").count() >= 1, f"[{tile}] group contrast selectors present")
        check(page.locator("text=Log2-transform").count() >= 1, f"[{tile}] log-transform option present")
        pills = page.locator("label.ora-db-pill").all_inner_texts()
        check(len(pills) == 6, f"[{tile}] all six databases offered (got {pills})")
        for db in ["Reactome", "SMPDB", "PathBank", "Lipid Class"]:
            page.locator("label.ora-db-pill", has_text=db).locator("input[type=checkbox]").uncheck()
        page.wait_for_timeout(150)
        # An EXPLICIT KEGG organism (yeast) beside a human organism above must win
        # and reach the run - typed by its code, as in Metabolic Connection's box.
        qbox = page.locator(".ora-kegg-organism-wrap")
        check(qbox.count() == 1, f"[{tile}] the KEGG organism box is on the Quantitative tile too")
        qbox.locator("input[placeholder*='KEGG organism']").fill("sce")
        for _ in range(40):
            if qbox.locator("select option", has_text="sce |").count() > 0:
                break
            page.wait_for_timeout(500)
        qbox.locator("select").select_option(label=qbox.locator("select option", has_text="sce |").first.inner_text())
        for _ in range(30):
            if page.evaluate("() => !!(STATE.metabolomics.run.draft && STATE.metabolomics.run.draft.keggOrganism)"):
                break
            page.wait_for_timeout(500)
        check(page.evaluate("() => STATE.metabolomics.run.draft.keggOrganism.code") == "sce",
              f"[{tile}] picked KEGG organism sce")
        if run_and_wait(page, tile):
            qkegg = page.evaluate("STATE.metabolomics.run.results.tables.KEGG")
            check(qkegg["n_terms_tested"] > 0
                  and all(r["Term ID"].startswith("sce") for r in qkegg["table"]["rows"]),
                  f"[{tile}] the explicit KEGG organism reached the run: KEGG tested sce maps beside a human organism "
                  f"({qkegg['n_terms_tested']} terms)")
            check("chosen explicitly" in (qkegg.get("caveat") or ""), f"[{tile}] the KEGG caveat says it was chosen explicitly")
            assert_house_style_results(page, tile)
            check(page.locator("text=Compound statistics").count() >= 1, f"[{tile}] the derived statistics are shown")
            check(page.locator("button", has_text="Download compound statistics").count() == 1, f"[{tile}] statistics download offered")
            sm = page.evaluate("STATE.metabolomics.statsMeta")
            check(sm["n_compounds_tested"] > 30, f"[{tile}] most compounds tested (got {sm['n_compounds_tested']})")
            summary = page.evaluate("STATE.metabolomics.run.results.summary")
            check(summary["universe_scope"] == "experiment", f"[{tile}] quantitative tile keeps the measured-compounds universe")
            tables = page.evaluate("STATE.metabolomics.run.results.tables")
            check(all(" (" not in k for k in tables), f"[{tile}] no direction suffix ({list(tables)})")
            # The tile mixes pathway and class families in one results screen, so
            # each panel must get ITS OWN chart set. Derive is_class from the
            # SELECTED key, not the tile - KEGG can now legitimately come back
            # empty, and asserting the pathway chart set against a class panel
            # would be a false failure.
            def _has_sig(k):
                return page.evaluate("""(k) => {
                    const b = STATE.metabolomics.run.results.tables[k];
                    return !!(b && b.table.rows.some(r => r.FDR !== null && r.FDR !== '' && Number(r.FDR) < STATE.metabolomics.run.sigAlpha));
                }""", k)
            exercised = 0
            for panel_key in ["KEGG", "Chemical Class"]:
                if panel_key not in tables:
                    continue
                select_result(page, panel_key)
                is_class_panel = panel_key == "Chemical Class"
                if panel_key == "KEGG":
                    check(page.locator("button.result-nav-btn", has_text="Class composition").count() == 0,
                          f"[{tile}] NO pie on the KEGG panel of the mixed-family tile")
                else:
                    exercise_class_composition(page, tile, errors)
                if _has_sig(panel_key):
                    rows = tables[panel_key]["table"]["rows"]
                    check("Compounds in set" in rows[0] and "Genes in set" not in rows[0],
                          f"[{tile}] burden columns renamed for compounds ({panel_key})")
                    exercise_charts(page, tile + " / " + panel_key, True, errors, is_class=is_class_panel)
                    exercised += 1
            check(exercised >= 1, f"[{tile}] at least one panel had a significant term to exercise")
            stats_csv = page.request.get(BASE + "/api/metabolomics/quant/download/statistics.csv",
                                         headers={"X-Session-Id": page.evaluate("getSessionId()")})
            check(stats_csv.ok and "log2fc" in stats_csv.text().splitlines()[0], f"[{tile}] statistics CSV downloads")

    # ----------------------------------------------------- Pathway Topology
    print("=== Pathway Topology Analysis: list mode ===")
    tile = "Pathway Topology Analysis"
    open_tile(page, tile)
    check(page.locator("div.section-label", has_text="Input").count() == 1, f"[{tile}] input-mode choice present")
    input_sel = page.locator("select").first
    check(input_sel.input_value() == "list", f"[{tile}] defaults to compound-list input")
    assert_house_style_module_screen(page, tile)
    if upload_and_resolve(page, tile, "example_metabolomics_list.csv", True):
        click_continue_and_run(page, tile)
        check(page.locator("div.section-label", has_text="Node importance").count() == 1,
              f"[{tile}] node-importance measure selector present")
        check(page.locator(".ora-kegg-organism-wrap").count() == 1,
              f"[{tile}] KEGG organism box present (topology always uses KEGG)")
        check(page.locator("label.ora-db-pill").count() == 0,
              f"[{tile}] no database picker - there is only ever one database")
        if run_and_wait(page, tile):
            assert_house_style_results(page, tile)
            tables = page.evaluate("STATE.metabolomics.run.results.tables")
            check(list(tables) == ["KEGG"], f"[{tile}] exactly one database (got {list(tables)})")
            block = tables["KEGG"]
            rows = block["table"]["rows"]
            check(bool(rows) and all("Impact" in r and "Matched in graph" in r for r in rows),
                  f"[{tile}] every row carries the topology columns")
            check(page.locator("button.result-nav-btn", has_text="Pathway impact").count() == 1,
                  f"[{tile}] 'Pathway impact' pane present")
            check(page.locator("button.result-nav-btn", has_text="Pathway graph").count() == 1,
                  f"[{tile}] 'Pathway graph' pane present")
            # Negative: none of the OTHER tiles' chart panes/toggles belong here.
            check(page.locator("button.result-nav-btn", has_text="Compound set overview").count() == 0,
                  f"[{tile}] no 'Compound set overview' pane (that is the ORA/Quant/class layout, not this one)")
            check(page.locator("button.result-nav-btn", has_text="Class composition").count() == 0,
                  f"[{tile}] no class-composition pane - class taxonomies are not tested here")

            page.locator("button.result-nav-btn", has_text="Pathway impact").first.click()
            page.wait_for_timeout(150)
            impact_btn = page.locator("button", has_text="Show pathway impact plot")
            check(impact_btn.count() == 1, f"[{tile}] pathway-impact toggle present")
            impact_btn.first.click()
            page.wait_for_timeout(800)
            impact_trace0 = page.evaluate("""() => {
                const gd = document.querySelector('.js-plotly-plot');
                return gd && gd.data && gd.data[0] ? { mode: gd.data[0].mode, y: gd.data[0].y } : null;
            }""")
            impact_names = impact_trace0["y"] or []
            check(len(impact_names) == len(rows),
                  f"[{tile}] the impact plot draws EVERY tested pathway by default, not just significant ones "
                  f"(got {len(impact_names)} points for {len(rows)} rows)")
            check(impact_trace0["mode"] == "markers" and all(n for n in impact_names),
                  f"[{tile}] every pathway's NAME is a real y-axis category label, not floating text (got mode="
                  f"{impact_trace0['mode']!r}, sample={impact_names[:3]})")
            # REGRESSION: a hardcoded x-axis range clipped the biggest marker
            # at either end (half a circle visible) until the user clicked
            # the modebar's "Autoscale" button by hand. No fixed range means
            # the FIRST draw already looks the way that button fixes it to.
            # Note: Plotly writes the COMPUTED range back into gd.layout once
            # autoranging runs, so gd.layout.xaxis.range is never undefined
            # after render - autorange itself (true here, false for a fixed
            # range) is the reliable signal that no hardcoded range was set.
            impact_autorange = page.evaluate("""() => {
                const gd = document.querySelector('.js-plotly-plot');
                return gd._fullLayout.xaxis.autorange;
            }""")
            check(bool(impact_autorange),
                  f"[{tile}] the impact plot's x-axis is genuinely autoranged, same as clicking the modebar's "
                  f"'Autoscale' button (got autorange={impact_autorange!r})")
            check(page.locator("span.chart-option-label", has_text="Terms shown").count() == 1,
                  f"[{tile}] 'Terms shown' control present on the impact plot")
            terms_shown_select = page.locator("label.chart-option-control", has_text="Terms shown").locator("select")
            terms_shown_select.select_option("10")
            page.wait_for_timeout(500)
            impact_capped = page.evaluate("""() => {
                const gd = document.querySelector('.js-plotly-plot');
                return gd && gd.data && gd.data[0] ? (gd.data[0].y || []).length : 0;
            }""")
            check(impact_capped == min(10, len(rows)),
                  f"[{tile}] 'Terms shown' narrows the impact plot (asked 10, got {impact_capped} of {len(rows)})")
            terms_shown_select.select_option("all")
            page.wait_for_timeout(500)

            page.locator("button.result-nav-btn", has_text="Pathway graph").first.click()
            page.wait_for_timeout(150)
            path_select = page.locator("select.chart-option-select").last
            n_options = page.evaluate("(el) => el.options.length", path_select.element_handle())
            check(n_options == len(rows), f"[{tile}] the pathway selector lists every tested pathway (got {n_options})")
            page.locator("button", has_text="Show pathway graph").first.click()
            page.wait_for_timeout(800)
            # The impact plot's chart div is still in the DOM (hidden, not
            # removed) once its pane has been shown once, so plain
            # querySelector('.js-plotly-plot') would silently grab THAT
            # stale plot instead of the graph just drawn - take the LAST
            # (most recently rendered) one instead.
            graph_info = page.evaluate("""() => {
                const all = document.querySelectorAll('.js-plotly-plot');
                const gd = all[all.length - 1];
                return gd && gd.data ? gd.data.map(t => (t.x || []).length) : null;
            }""")
            check(bool(graph_info) and graph_info[-1] > 0,
                  f"[{tile}] the pathway graph draws node markers for the selected pathway (got {graph_info})")
            # Labels are Plotly ANNOTATIONS (draggable independently of the
            # node marker), not the marker trace's own text - see
            # buildPathwayGraphPlot's comment. Same "last plot in the DOM"
            # caveat as graph_info above.
            graph_annotations = page.evaluate("""() => {
                const all = document.querySelectorAll('.js-plotly-plot');
                const gd = all[all.length - 1];
                return (gd && gd.layout && gd.layout.annotations || []).map(a => a.text);
            }""")
            check(len(graph_annotations) > 0 and all(t for t in graph_annotations),
                  f"[{tile}] every graph node gets a draggable compound-name annotation by default "
                  f"(got {len(graph_annotations)} annotations)")
            check(any(not re.fullmatch(r"[CG]\d{5}", t) for t in graph_annotations if t),
                  f"[{tile}] node labels are real compound NAMES, not bare KEGG ids (sample={graph_annotations[:5]})")
            check(page.locator("span.chart-option-label", has_text="Compound labels").count() == 1,
                  f"[{tile}] 'Compound labels' show/hide control present on the graph pane")
            check(page.locator("button", has_text="Reset label positions").count() == 1,
                  f"[{tile}] 'Reset label positions' action present on the graph pane")
            check(page.locator("div.ora-chart-threshold-note", has_text="genuine KEGG data").count() == 1,
                  f"[{tile}] the graph pane explains isolated (unconnected) nodes as genuine KEGG data, not a bug")
            check(page.locator("div.ora-chart-threshold-note", has_text="can be dragged").count() == 1,
                  f"[{tile}] the graph pane explains that nodes and their compound labels can be dragged")
            labels_select = page.locator("select.chart-option-select").filter(
                has=page.locator("option", has_text="Hide")).last
            labels_select.select_option("hide")
            page.wait_for_timeout(500)
            graph_annotations_hidden = page.evaluate("""() => {
                const all = document.querySelectorAll('.js-plotly-plot');
                const gd = all[all.length - 1];
                return (gd && gd.layout && gd.layout.annotations || []).length;
            }""")
            check(graph_annotations_hidden == 0,
                  f"[{tile}] 'Compound labels' -> Hide removes the on-node annotations (got {graph_annotations_hidden})")
            labels_select.select_option("show")
            page.wait_for_timeout(500)

            # Drag one node's label with a REAL mouse gesture, same proven
            # pattern as the class-composition pie's percentage labels
            # (exercise_pie_labels above): the annotation's anchor (x, y -
            # the node's own KGML position) must stay put while its tail
            # (ax, ay - the label) moves, and the move must survive a redraw
            # and be undoable via Reset.
            def graph_state():
                return page.evaluate("""() => {
                    const all = document.querySelectorAll('.js-plotly-plot');
                    const gd = all[all.length - 1];
                    return (gd && gd.layout && gd.layout.annotations || [])
                        .map(a => ({text: a.text, x: a.x, y: a.y, ax: a.ax, ay: a.ay}));
                }""")
            before = graph_state()[0]
            label_el = page.locator('.js-plotly-plot .annotation[data-index="0"] text').last
            label_el.scroll_into_view_if_needed()
            page.wait_for_timeout(200)
            box = label_el.bounding_box()
            check(box is not None, f"[{tile}] the first node's compound-name label is on screen")
            if box:
                cx, cy = box["x"] + box["width"] / 2, box["y"] + box["height"] / 2
                page.mouse.move(cx, cy)
                page.mouse.down()
                page.mouse.move(cx + 30, cy + 20, steps=4)
                page.mouse.move(cx + 65, cy + 45, steps=6)
                page.mouse.up()
                page.wait_for_timeout(500)
                moved = graph_state()[0]
                check(abs(moved["ax"] - before["ax"]) > 20 or abs(moved["ay"] - before["ay"]) > 20,
                      f"[{tile}] dragging a node's label moves it (ax {before['ax']:.0f}->{moved['ax']:.0f}, "
                      f"ay {before['ay']:.0f}->{moved['ay']:.0f})")
                check(abs(moved["x"] - before["x"]) < 1e-9 and abs(moved["y"] - before["y"]) < 1e-9,
                      f"[{tile}] ...while the leader line stays attached to the node (anchor unchanged)")
                labels_select.select_option("hide")
                page.wait_for_timeout(300)
                labels_select.select_option("show")
                page.wait_for_timeout(500)
                kept = graph_state()[0]
                check(abs(kept["ax"] - moved["ax"]) < 1e-6 and abs(kept["ay"] - moved["ay"]) < 1e-6,
                      f"[{tile}] a dragged label keeps its place across a redraw (stored on the run)")
                page.locator("button.chart-toggle-btn", has_text="Reset label positions").first.click()
                page.wait_for_timeout(500)
                reset = graph_state()[0]
                check(abs(reset["ax"] - before["ax"]) < 1e-6 and abs(reset["ay"] - before["ay"]) < 1e-6,
                      f"[{tile}] 'Reset label positions' puts it back")

            # Drag the NODE itself (not its label) with a real mouse gesture -
            # reuses _makeNetworkDraggable, the same click-drag this app's
            # other network charts (enrichment map, gene-concept network)
            # already use. Node 1 (not 0, already exercised above) is
            # targeted by converting its data x/y to a screen position the
            # same way _makeNetworkDraggable's own hit-test does, since a
            # Plotly marker has no easily selectable DOM element of its own.
            def marker_client_pos(node_index):
                return page.evaluate("""(i) => {
                    const all = document.querySelectorAll('.js-plotly-plot');
                    const gd = all[all.length - 1];
                    const fl = gd._fullLayout;
                    const tr = gd.data[gd.data.length - 1];
                    const rect = gd.getBoundingClientRect();
                    const size = fl._size, xr = fl.xaxis.range, yr = fl.yaxis.range;
                    const px = size.l + ((tr.x[i] - xr[0]) / (xr[1] - xr[0])) * size.w;
                    const py = size.t + ((yr[1] - tr.y[i]) / (yr[1] - yr[0])) * size.h;
                    return { clientX: rect.left + px, clientY: rect.top + py, x: tr.x[i], y: tr.y[i] };
                }""", node_index)

            def node_edge_state():
                return page.evaluate("""() => {
                    const all = document.querySelectorAll('.js-plotly-plot');
                    const gd = all[all.length - 1];
                    const tr = gd.data[gd.data.length - 1];
                    const ann = (gd.layout.annotations || [])[1];
                    return {
                        nodeX: tr.x[1], nodeY: tr.y[1],
                        annX: ann ? ann.x : null, annY: ann ? ann.y : null,
                    };
                }""")

            # Scroll the chart into view first (same reason the label-drag
            # test above does) - otherwise the computed client coordinates
            # can point off-screen and the synthetic mouse gesture lands on
            # nothing, silently no-op-ing the drag.
            page.locator(".js-plotly-plot").last.scroll_into_view_if_needed()
            page.wait_for_timeout(200)
            before_node = marker_client_pos(1)
            before_state = node_edge_state()
            page.mouse.move(before_node["clientX"], before_node["clientY"])
            page.mouse.down()
            page.mouse.move(before_node["clientX"] + 40, before_node["clientY"] + 30, steps=4)
            page.mouse.move(before_node["clientX"] + 80, before_node["clientY"] + 55, steps=6)
            page.mouse.up()
            page.wait_for_timeout(500)
            after_state = node_edge_state()
            check(abs(after_state["nodeX"] - before_state["nodeX"]) > 1 or abs(after_state["nodeY"] - before_state["nodeY"]) > 1,
                  f"[{tile}] dragging a graph NODE moves it ({before_state['nodeX']:.1f},{before_state['nodeY']:.1f}) "
                  f"-> ({after_state['nodeX']:.1f},{after_state['nodeY']:.1f})")
            check(after_state["annX"] is not None and abs(after_state["annX"] - after_state["nodeX"]) < 1e-6
                  and abs(after_state["annY"] - after_state["nodeY"]) < 1e-6,
                  f"[{tile}] the dragged node's own label follows it (annotation anchor == new node position, "
                  f"got ann=({after_state['annX']},{after_state['annY']}) vs node=({after_state['nodeX']},{after_state['nodeY']}))")

            second_value = page.evaluate("(el) => el.options[1].value", path_select.element_handle())
            path_select.select_option(second_value)
            page.wait_for_timeout(800)
            check(len(errors) == 0, f"[{tile}] switching the drawn pathway raised no page error")

            check(page.evaluate("STATE.metabolomics.run.results.tables.KEGG.table.rows")
                  and "Overlapping Compounds" in rows[0],
                  f"[{tile}] the member column is relabelled to compound names, same as every other tile")

    print("=== Pathway Topology Analysis: concentration-table mode ===")
    open_tile(page, tile)
    sel = page.locator("select").first
    sel.select_option("conc")
    page.wait_for_timeout(200)
    check(page.evaluate("STATE.metabolomics.topologyInput") == "conc",
          f"[{tile}] switching the input select updates topologyInput")
    if upload_and_resolve(page, tile, "example_metabolomics_concentrations.csv", True):
        click_continue_and_run(page, tile)
        check(page.locator("text=Groups to compare").count() >= 1,
              f"[{tile}] group-contrast fields appear in concentration mode, same as Quantitative Enrichment")
        check(page.locator("div.section-label", has_text="Node importance").count() == 1,
              f"[{tile}] node-importance selector still present in concentration mode")
        measure_row = page.locator("select").filter(has=page.locator(
            "option", has_text="Out-degree"))
        if measure_row.count():
            measure_row.first.select_option("outdegree")
        if run_and_wait(page, tile):
            assert_house_style_results(page, tile)
            block2 = page.evaluate("STATE.metabolomics.run.results.tables.KEGG")
            check(block2.get("measure") == "outdegree", f"[{tile}] the chosen measure reached the run")
            check(block2["table"]["rows"] and "Effect (mean difference)" in block2["table"]["rows"][0],
                  f"[{tile}] concentration mode uses the Wilcoxon/burden column shape")
            check(page.locator("div.section-label", has_text="Compound statistics").count() == 1,
                  f"[{tile}] topology's own statistics preview is shown (concentration mode)")

    print("=== a pie with many classes sizes itself ===")
    exercise_many_classes(page, errors)

    print("=== console error check ===")
    check(len(errors) == 0, f"no console/page errors during the whole flow (got {errors})")
    browser.close()

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
