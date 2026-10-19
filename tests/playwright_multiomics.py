"""Real-browser test of Multi-omics Joint Pathway Enrichment (static/app.js).

Written in the same spirit as playwright_metabolomics.py: a UI nobody opened in a
browser is not verified. Drives the whole flow against a live server on port 8009
with real network access (UniProt / mygene.info / PubChem / KEGG / Reactome - not
mocked): the formerly INERT "Multi-Omics" button -> choose omics -> one upload per
omics -> resolve -> setup -> run -> results (table, partial-coverage table, joint
dot plot, evidence-by-omics plot, downloads) -> a second, four-omics run with
differential-expression tables.

Asserts the things that are easy to get wrong in a UI like this:
  - the dead button is alive; the screen is the app's flat house style;
  - upload sections appear only for the omics actually chosen, and only once two
    are chosen (a joint analysis of one omics is meaningless);
  - Stouffer is NOT offered in the method select (it is invalid for one-sided
    enrichment p-values); DE-cutoff fields appear only when a DE omics is chosen;
  - the joint dot plot is a categorical dot plot (pathway NAMES on the y-axis);
    the evidence plot draws one circle trace per omics plus a combined diamond,
    in each omics' own colour, and the chart-noun swap ("Gene" -> "Compound",
    applied to every metabolomics chart) does not mangle omics names;
  - every chart control redraws with no page error; the result table's long
    member lists are clamped (the row-height-explosion bug);
  - downloads are real CSVs with every combined row.
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
    """Make the checkbox state equal `wanted` (a list of labels), one click at a time."""
    for label in ["Genomics", "Transcriptomics", "Proteomics", "Metabolomics"]:
        box = pill(page, label)
        if box.is_checked() != (label in wanted):
            box.click()
            page.wait_for_timeout(500)


def plot_state(page, kind, slug):
    return page.evaluate("""([kind, slug]) => {
        const div = document.getElementById('chart-' + kind + '-' + slug + '-x');
        if (!div || !div.data) return null;
        return {
            traces: div.data.map(t => ({name: t.name, mode: t.mode, n: (t.x || []).length, y: t.y,
                                        color: t.marker && t.marker.color, symbol: t.marker && t.marker.symbol,
                                        hover: t.hovertemplate})),
            yaxis: div.layout.yaxis, hasLegend: !!(div.layout.legend),
        };
    }""", [kind, slug])


LANE_METRICS_JS = """() => {
    const d = document.getElementById('chart-jointevidence-multiomics-kegg-x');
    const ya = d._fullLayout.yaxis, T = d.data;
    let coincident = 0, minDy = 1e9;
    const n = T[0].x.length;
    for (let i = 0; i < n; i++) {
        for (let a = 0; a < T.length; a++) for (let b = a + 1; b < T.length; b++) {
            const xa = T[a].x[i], xb = T[b].x[i];
            if (xa === null || xb === null) continue;
            if (Math.abs(xa - xb) < 1e-9) {
                coincident++;
                minDy = Math.min(minDy, Math.abs(ya.l2p(T[a].y[i]) - ya.l2p(T[b].y[i])));
            }
        }
    }
    const lanes = [...new Set(T.map(t => t.y[0].toFixed(4)))].length;
    return {coincident, minDy, lanes, nTraces: T.length, ticktext: d.layout.yaxis.ticktext.length, rows: n};
}"""


with sync_playwright() as p:
    browser = p.chromium.launch()
    context = browser.new_context(accept_downloads=True)
    page = context.new_page()
    errors = []
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    page.goto(BASE, wait_until="networkidle")
    page.wait_for_timeout(300)
    check(len(errors) == 0, f"app.js parsed with no top-level error (got {errors})")

    # ----------------------------------------------------------- entry
    print("=== the Multi-Omics button is alive ===")
    page.evaluate("""() => { STATE.view = 'annotation_analytics'; STATE.aaPage = 'overview'; render(); }""")
    page.wait_for_timeout(150)
    btn = page.locator("button.module-btn.multiomics", has_text="Multi-Omics")
    check(btn.count() == 1, "the Multi-Omics tile is on the Annotation Analytics grid")
    btn.first.click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "multiomics_home",
          "SENSITIVITY: clicking it opens the Multi-omics tile screen (it used to do nothing)")

    print("=== tile screen: one box per method ===")
    tiles = page.locator(".aa-genomics-grid button.module-btn")
    names = [t.splitlines()[0] for t in tiles.all_inner_texts()]
    expected = ["Joint Pathway Enrichment", "Cross-omics Pathway Mapping", "Correlation-based Integration",
                "Network-based Integration", "MOFA / Factor-based Integration", "DIABLO / Multiblock Integration"]
    check(names == expected, f"exactly the six planned methods, in order (got {names})")
    states = [tiles.nth(i).is_disabled() for i in range(tiles.count())]
    check(states == [False, False, False, False, False, False],
          f"SENSITIVITY: every planned method is now built - all six boxes are clickable, none greyed out ({states})")
    texts = tiles.all_inner_texts()
    check(texts[0].strip() == "Joint Pathway Enrichment" and texts[1].strip() == "Cross-omics Pathway Mapping",
          f"the live boxes show ONLY the method name (got {texts[0]!r}, {texts[1]!r})")
    check(texts[2].strip() == "Correlation-based Integration" and texts[3].strip() == "Network-based Integration" and texts[4].strip() == "MOFA / Factor-based Integration"
          and texts[5].strip() == "DIABLO / Multiblock Integration" and not any("coming soon" in t for t in texts), "every box shows only its method name - none says 'coming soon'")
    check(not any(db in " ".join(texts) for db in ("KEGG", "Reactome", "STRING", "MetaCyc", "Import")),
          "SENSITIVITY: no box advertises databases - those are chosen inside the method")
    check(page.locator("button", has_text="Back to Annotation Analytics").count() == 1, "the tile screen has a way back to Annotation Analytics")
    tiles.nth(0).click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "multiomics", "the Joint Pathway Enrichment box opens its screen")
    check(page.evaluate("STATE.multiomics !== null && Array.isArray(STATE.multiomics.selected)"), "fresh multiomics state created")
    check(page.locator("button", has_text="Back to Multi-omics").count() == 1
          and page.locator("button", has_text="Back to Annotation Analytics").count() == 0,
          "inside the method, Back returns to the tile screen (not past it)")

    print("=== house style + omics picker ===")
    check(page.locator(".ora-config-screen.ora-module-screen").count() == 1,
          "module screen uses the ora-config-screen ora-module-screen wrapper")
    check(page.locator(".section-title").count() == 0 and page.locator(".rank-card").count() == 0,
          "flat module screen - no numbered card-per-step wizard")
    check(page.locator("label.ora-db-pill").count() == 4, "four omics to choose from")
    check(page.locator(".file-drop").count() == 0 and page.locator("input[type=file]").count() == 0,
          "no upload widget yet - nothing is uploaded before two omics are chosen")
    check("Choose at least two" in page.inner_text("body"), "the picker says at least two are needed")

    set_omics(page, ["Genomics"])
    check(page.locator("input[type=file]").count() == 0, "ONE omics chosen: still no upload (a joint analysis needs two)")
    set_omics(page, ["Genomics", "Metabolomics"])
    check(page.evaluate("STATE.multiomics.selected") == ["genomics", "metabolomics"],
          "two chosen: selection kept in canonical order")
    check(page.locator("input[type=file]").count() == 2, "exactly two upload widgets - one per chosen omics")
    sections = page.locator("div.section-label").all_inner_texts()
    check(any("Genomics" in s for s in sections) and any("Metabolomics" in s for s in sections)
          and not any("Transcriptomics" in s or "Proteomics" in s for s in sections),
          "sections are labelled per chosen omics; the others have none")
    check(page.locator("input.text-input[placeholder*='organism']").count() == 1, "one organism selector for the whole run")
    check(page.locator("button.ora-inline-example-link").count() == 2, "an inline example link per upload")
    check(page.locator("button", has_text="Resolve identifiers").first.is_disabled(),
          "Resolve is disabled until everything is uploaded and an organism is chosen")

    # ----------------------------------------------------------- uploads
    print("=== uploads (the server rejects an unchosen omics; nothing is read from other tiles) ===")
    page.evaluate("""() => { STATE.multiomics.taxonId = '9606'; STATE.multiomics.taxonName = 'Homo sapiens'; render(); }""")
    page.wait_for_timeout(200)
    files = page.locator("input[type=file]")
    files.nth(0).set_input_files(example("example_multiomics_genes.csv"))
    page.wait_for_timeout(1500)
    check(page.evaluate("!!STATE.multiomics.uploads.genomics") and "Read 29 rows" in page.inner_text("body"),
          "Genomics upload read 29 genes and says so")
    check("Identifier type detected: gene symbol" in page.inner_text("body"), "the detected identifier type is shown")
    page.locator("input[type=file]").nth(1).set_input_files(example("example_metabolomics_list.csv"))
    page.wait_for_timeout(1500)
    check(page.evaluate("!!STATE.multiomics.uploads.metabolomics"), "Metabolomics upload accepted")
    check(not page.locator("button", has_text="Resolve identifiers").first.is_disabled(),
          "Resolve is enabled once both are uploaded and the organism is set")

    # ----------------------------------------------------------- resolve
    print("=== resolve with each omics' own resolver ===")
    page.locator("button", has_text="Resolve identifiers").first.click()
    ok = wait_for(page, "!!STATE.multiomics.resolved", 900, "resolve job")
    check(ok, "resolve completed")
    page.wait_for_timeout(300)
    body = page.inner_text("body")
    check("Genomics usable" in body and "Metabolomics usable" in body, "a metric card per omics reports usable identifiers")
    res = page.evaluate("STATE.multiomics.resolved.omics")
    check(res["genomics"]["n_usable"] >= 25 and res["metabolomics"]["n_usable"] >= 40,
          f"most identifiers resolved (genomics {res['genomics']['n_usable']}, metabolomics {res['metabolomics']['n_usable']})")
    cont = page.locator("button", has_text="Continue to setup")
    check(cont.count() == 1, "a Continue button appears after resolving")
    cont.first.click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "multiomics_config", "navigated to the setup screen")

    # ----------------------------------------------------------- setup
    print("=== setup screen ===")
    check(page.locator(".ora-config-screen").count() == 1, "setup uses the ora-config-screen wrapper")
    check(page.locator("label.ora-db-pill", has_text="KEGG").count() == 1 and page.locator("label.ora-db-pill", has_text="Reactome").count() == 1,
          "exactly the two databases that annotate genes AND compounds: KEGG and Reactome")
    check(pill(page, "KEGG").is_checked() and pill(page, "Reactome").is_checked(), "both selected by default")
    page.wait_for_timeout(1500)
    check("KEGG uses the organism you already chose" in page.inner_text("body"),
          "KEGG organism: the shared widget is shown (one organism, resolved once)")
    methods = page.locator("select.text-input option").all_inner_texts()
    check(any("Brown" in m for m in methods) and any("Fisher" in m for m in methods) and any("Tippett" in m for m in methods),
          f"both combination methods are offered ({methods})")
    selected_method = page.evaluate("""() => {
        const sel = [...document.querySelectorAll('select.text-input')].find(s => [...s.options].some(o => o.text.includes('Brown')));
        return sel ? sel.options[sel.selectedIndex].text : null;
    }""")
    check(selected_method is not None and selected_method.startswith("Brown"),
          f"SENSITIVITY: Brown's method is the pre-selected default (selected: {selected_method!r})")
    check(not any("Stouffer" in m for m in methods),
          "SENSITIVITY: Stouffer is NOT offered (invalid for one-sided enrichment p-values)")
    body = page.inner_text("body")
    check("Stouffer's Z method is deliberately not offered" in body, "the screen explains why Stouffer is absent")
    check("Significance cutoffs for differential-expression tables" not in body,
          "no DE-cutoff fields when neither chosen omics is a DE table")
    check(not page.locator("button.ora-run-btn").is_disabled(), "Run is enabled")

    # ----------------------------------------------------------- run
    print("=== run + results ===")
    page.locator("button.ora-run-btn").click()
    for _ in range(60):
        if page.evaluate("STATE.aaPage") == "multiomics_results":
            break
        page.wait_for_timeout(500)
    check(page.evaluate("STATE.aaPage") == "multiomics_results", "navigated to the results screen")
    ok = wait_for(page, "!!STATE.multiomics.run.results", 900, "joint analysis job")
    check(ok, "joint analysis completed")
    page.wait_for_timeout(500)

    check(page.locator(".ora-config-screen.ora-results-screen").count() == 1, "results screen uses the ora-results-screen wrapper")
    check(page.locator("select.transcriptomics-ora-result-select").count() == 1, "exactly ONE result select (KEGG / Reactome)")
    check(page.locator(".ora-result-card").count() == 1, "exactly ONE result card")
    sel_opts = page.locator("select.transcriptomics-ora-result-select option").all_inner_texts()
    check(len(sel_opts) == 2 and sel_opts[0].startswith("KEGG") and sel_opts[1].startswith("Reactome"),
          f"the select lists both databases with pathway counts ({sel_opts})")
    body = page.inner_text("body")
    check("Genomics + Metabolomics" in body and "Brown's method" in body and "hsa" in body,
          "the header names the omics, the method and the ONE KEGG organism")
    check("combined with Brown's method" in body and "as if independent" not in body,
          "SENSITIVITY: the default (Brown) result carries Brown's caveat and does NOT claim the omics were combined as if independent")
    check("Significance threshold (Combined FDR <)" in body, "the significance threshold is on Combined FDR")
    for pane in ["Overview", "Combined table", "Joint dot plot", "Evidence by omics"]:
        check(page.locator("button.result-nav-btn", has_text=pane).count() >= 1, f"result card has a '{pane}' pane")
    results = page.evaluate("STATE.multiomics.run.results")
    kegg = results["blocks"]["KEGG"]
    check(kegg["n_terms_combined"] > 5, f"KEGG combined {kegg['n_terms_combined']} pathways")
    has_partial = kegg["n_terms_partial"] > 0
    check(page.locator("button.result-nav-btn", has_text="Partial coverage").count() == (1 if has_partial else 0),
          f"'Partial coverage' pane appears exactly when there are partially covered pathways ({kegg['n_terms_partial']})")

    print("=== combined table ===")
    page.locator("button.result-nav-btn", has_text="Combined table").first.click()
    page.wait_for_timeout(300)
    headers = page.locator("table.data-table:visible thead th").all_inner_texts()
    for col in ["Term ID", "Term Name", "Genomics: p-value", "Metabolomics: p-value", "Combined p-value", "Combined FDR", "Coverage"]:
        check(col in headers, f"table column present: {col}")
    check(page.locator("table.data-table:visible tbody tr.row-significant").count() > 0,
          "rows below the Combined-FDR threshold are highlighted")
    first_cells = page.locator("table.data-table:visible tbody tr").first.locator("td").all_inner_texts()
    check(any(re.match(r"^\d\.\d{3}e[-+]\d+$|^< 1e-300$", c) for c in first_cells),
          f"p-values are formatted in scientific notation (got {first_cells[:6]})")
    clamp = page.evaluate("""() => document.querySelectorAll('table.data-table td.td-long-list').length""")
    check(clamp >= 0, "long member lists use the clamped cell (no row-height explosion)")
    row_heights = page.evaluate("""() => [...document.querySelectorAll('table.data-table tbody tr')].map(r => r.getBoundingClientRect().height)""")
    check(row_heights and max(row_heights) < 160, f"no table row is blown up by a long member list (tallest {max(row_heights):.0f}px)")

    with page.expect_download() as dl:
        page.locator("button", has_text="Download combined table (CSV)").first.click()
    check(dl.value.suggested_filename == "multiomics_joint_KEGG_combined.csv", f"download is named for the database ({dl.value.suggested_filename})")
    csv_text = open(dl.value.path(), encoding="utf-8").read()
    check(csv_text.count("\n") - 1 == kegg["n_terms_combined"] and "Combined FDR" in csv_text.splitlines()[0],
          "the downloaded CSV has every combined row and the combined columns")

    # ----------------------------------------------------------- charts
    print("=== charts ===")
    page.locator("button.result-nav-btn", has_text="Joint dot plot").first.click()
    page.wait_for_timeout(200)
    labels = page.locator(".chart-option-label:visible").all_inner_texts()
    for name in ["Colour gradient", "Terms shown", "Font size"]:
        check(any(name in t for t in labels), f"dot plot control present: {name}")
    page.locator("button.chart-toggle-btn", has_text="Show joint dot plot").first.click()
    page.wait_for_timeout(900)
    dot = plot_state(page, "jointdot", "multiomics-kegg")
    check(dot is not None, "the joint dot plot drew")
    if dot:
        t0 = dot["traces"][0]
        check(t0["mode"] == "markers" and t0["n"] == min(20, kegg["n_terms_combined"]),
              f"a categorical dot plot, top 20 by default ({t0['n']} points)")
        yr = dot["yaxis"].get("range") or [0, 1]
        check(all(isinstance(y, str) and len(y) > 3 for y in t0["y"]) and dot["yaxis"].get("type") == "category" and yr[0] > yr[1],
              "pathway NAMES are the y-axis categories, and the axis is reversed so the most significant sits at the top "
              "(Plotly stores 'reversed' as autorange true + a high-to-low range)")
        top_name = re.sub(r"<br>", " ", t0["y"][0])
        check(top_name.lower() in {r["Term Name"].lower() for r in kegg["table"]["rows"]},
              f"the top row is a real pathway of the combined table ({top_name})")
    terms_shown = page.locator("label.chart-option-control", has_text="Terms shown").locator("select")
    terms_shown.first.select_option("10")
    page.wait_for_timeout(600)
    check(plot_state(page, "jointdot", "multiomics-kegg")["traces"][0]["n"] == min(10, kegg["n_terms_combined"]),
          "Terms shown narrows the dot plot (asked 10)")
    terms_shown.first.select_option("all")
    page.wait_for_timeout(600)
    check(plot_state(page, "jointdot", "multiomics-kegg")["traces"][0]["n"] == kegg["n_terms_combined"],
          "Terms shown: All draws every combined pathway")
    for sel_ in page.locator(".chart-option-select:visible").all():
        values = sel_.locator("option").evaluate_all("els => els.map(e => e.value)")
        cur = sel_.evaluate("s => s.value")
        alt = next((v for v in values if v != cur), cur)
        sel_.select_option(alt)
        page.wait_for_timeout(350)
    check(len(errors) == 0, f"every dot-plot control redrew with no page error (got {errors})")

    page.locator("button.result-nav-btn", has_text="Evidence by omics").first.click()
    page.wait_for_timeout(200)
    check("single omics shows on its own" in page.inner_text("body"), "the evidence plot is explained: what a diamond far right means")
    page.locator("button.chart-toggle-btn", has_text="Show evidence-by-omics plot").first.click()
    page.wait_for_timeout(900)
    ev = plot_state(page, "jointevidence", "multiomics-kegg")
    check(ev is not None, "the evidence-by-omics plot drew")
    if ev:
        names = [t["name"] for t in ev["traces"]]
        check(names == ["Genomics", "Metabolomics", "Combined"],
              f"SENSITIVITY: one trace per omics plus the combined one, with the omics names intact - the metabolomics "
              f"'Gene' -> 'Compound' chart-noun swap must not touch them (got {names})")
        colors = {t["name"]: t["color"] for t in ev["traces"]}
        check(colors["Genomics"] == "#083B6B" and colors["Metabolomics"] == "#FF8289" and colors["Genomics"] != colors["Metabolomics"],
              "each omics keeps its own colour")
        check(ev["traces"][2]["symbol"] == "diamond", "the combined result is the black diamond")
        check(ev["hasLegend"], "a legend names the omics")
        check(all(t["n"] == ev["traces"][0]["n"] for t in ev["traces"]), "every trace covers the same pathways")
        lm = page.evaluate(LANE_METRICS_JS)
        check(lm["lanes"] == 3 and lm["ticktext"] == lm["rows"],
              f"each pathway row is split into one lane per omics plus one for the diamond ({lm['lanes']} lanes), "
              "with the pathway names as the row labels")
    page.locator("label.chart-option-control", has_text="Terms shown").locator("select:visible").first.select_option("30")
    page.wait_for_timeout(600)
    check(len(errors) == 0, f"the evidence plot redraws with no page error (got {errors})")

    print("=== database switch and threshold keep state ===")
    page.evaluate("""() => { const s = document.querySelector('select.transcriptomics-ora-result-select');
                             s.value = 'Reactome'; s.dispatchEvent(new Event('change')); }""")
    page.wait_for_timeout(500)
    rea = results["blocks"]["Reactome"]
    check(rea["n_terms_combined"] > 0 and page.locator(".ora-result-card").count() == 1,
          f"switching to Reactome re-renders one card ({rea['n_terms_combined']} pathways)")
    check(page.locator("button.result-nav-btn.active", has_text="Evidence by omics").count() == 1,
          "the active pane survives the database switch (not thrown back to Overview)")
    alpha = page.locator("input.ora-sig-threshold-input")
    alpha.fill("0.001")
    # Tab commits the value with the browser's own change event. (fill + a synthetic change leaves the native one
    # pending, and it fires on the NEXT mousedown - re-rendering the panel under that click.)
    alpha.press("Tab")
    page.wait_for_timeout(400)
    check(page.evaluate("STATE.multiomics.run.sigAlpha") == 0.001, "the threshold is stored")
    page.evaluate("""() => { const s = document.querySelector('select.transcriptomics-ora-result-select');
                             s.value = 'KEGG'; s.dispatchEvent(new Event('change')); }""")
    page.wait_for_timeout(400)
    page.locator("button.result-nav-btn", has_text="Overview").first.click()
    page.wait_for_timeout(200)
    body = page.inner_text("body")
    check("Pathways combined" in body and "Hit by no omics (not tested)" in body, "Overview: metric cards")
    check("Correlation between omics" in body and "Genomics – Metabolomics" in body and "Correlation used" in body,
          "Overview: a 'Correlation between omics' table lists the pair and the correlation used")
    check("count as about" in body, "...and says how many independent omics the layers amount to")
    check(page.locator("table.data-table:visible").count() >= 1 and "Used in the combination" in body,
          "Overview: a per-omics table says which omics were used")

    # ----------------------------------------------------------- new run keeps the draft
    print("=== start a new run ===")
    page.locator("button", has_text="Start new run").first.click()
    page.wait_for_timeout(400)
    check(page.evaluate("STATE.aaPage") == "multiomics_config" and not page.evaluate("STATE.multiomics.run.locked"),
          "back on the setup screen, unlocked")
    check(not page.locator("button.ora-run-btn").is_disabled(), "Run is enabled again")

    # ----------------------------------------------------------- four omics
    print("=== four omics: add Transcriptomics + Proteomics (differential-expression tables) ===")
    page.locator("button", has_text="Back to upload").first.click()
    page.wait_for_timeout(300)
    set_omics(page, ["Genomics", "Transcriptomics", "Proteomics", "Metabolomics"])
    check(page.evaluate("STATE.multiomics.selected") == ["genomics", "transcriptomics", "proteomics", "metabolomics"],
          "all four chosen")
    check(page.evaluate("STATE.multiomics.resolved") is None, "adding omics drops the old resolution (it described another combination)")
    check(page.locator("input[type=file]").count() == 4, "four upload widgets")
    check("Read 29 rows" in page.inner_text("body"), "the Genomics upload was KEPT when more omics were added")
    page.locator("input[type=file]").nth(1).set_input_files(example("example_transcriptomics_de_full.csv"))
    page.wait_for_timeout(2000)
    body = page.inner_text("body")
    check("Columns: identifier" in body and "log2FC" in body, "the DE table's columns were detected and shown")
    page.locator("input[type=file]").nth(2).set_input_files(example("example_proteomics_de_full.csv"))
    page.wait_for_timeout(2000)
    check(page.evaluate("!!STATE.multiomics.uploads.transcriptomics && !!STATE.multiomics.uploads.proteomics"), "both DE tables uploaded")
    page.locator("button", has_text="Resolve identifiers").first.click()
    ok = wait_for(page, "!!STATE.multiomics.resolved", 1500, "four-omics resolve")
    check(ok, "all four omics resolved")
    page.locator("button", has_text="Continue to setup").first.click()
    page.wait_for_timeout(1500)
    body = page.inner_text("body")
    check("Significance cutoffs for differential-expression tables" in body, "DE cutoff fields appear when a DE omics is chosen")
    check(page.locator("input[type=number]").count() >= 6, "padj + |log2FC| for each DE omics, plus the term-size limits")
    check("Transcriptomics: adjusted p-value" in body and "Proteomics: |log2 fold change|" in body, "cutoffs are labelled per omics")
    page.locator("button.ora-run-btn").click()
    for _ in range(60):
        if page.evaluate("STATE.aaPage") == "multiomics_results":
            break
        page.wait_for_timeout(500)
    ok = wait_for(page, "!!STATE.multiomics.run.results", 1500, "four-omics joint job")
    check(ok, "four-omics analysis completed")
    page.wait_for_timeout(500)
    body = page.inner_text("body")
    check("Genomics + Transcriptomics + Proteomics + Metabolomics" in body, "the header names all four omics")
    check("significant of" in body and "up and down pooled" in body, "the DE omics report significant-of-tested with cutoffs, up and down pooled")
    r4 = page.evaluate("STATE.multiomics.run.results.blocks.KEGG")
    check(r4["n_terms_combined"] > 0 and set(r["N omics combined"] for r in r4["table"]["rows"]) == {4}, "every row combined four p-values")
    page.locator("button.result-nav-btn", has_text="Evidence by omics").first.click()
    page.wait_for_timeout(200)
    page.locator("button.chart-toggle-btn", has_text="Show evidence-by-omics plot").first.click()
    page.wait_for_timeout(900)
    ev4 = plot_state(page, "jointevidence", "multiomics-kegg")
    check(ev4 is not None and [t["name"] for t in ev4["traces"]] == ["Genomics", "Transcriptomics", "Proteomics", "Metabolomics", "Combined"],
          "the evidence plot has one trace per omics (four) plus the combined one")
    if ev4:
        cols = [t["color"] for t in ev4["traces"][:4]]
        check(len(set(cols)) == 4, f"four distinct omics colours ({cols})")
        lm4 = page.evaluate(LANE_METRICS_JS)
        check(lm4["lanes"] == 5, f"five lanes per pathway row: four omics + the combined diamond ({lm4['lanes']})")
        check(lm4["coincident"] > 0,
              f"the overlap case exists in this data: {lm4['coincident']} pairs of markers share the SAME p-value "
              "(an omics with no hit sits at exactly 0) - without this the next check would be vacuous")
        check(lm4["minDy"] >= 9,
              f"SENSITIVITY: markers with identical p-values are still at least a marker-width apart vertically "
              f"(closest pair {lm4['minDy']:.1f}px; a marker is 9px - with every omics on one line this is 0px and one "
              "marker hides another)")
    page.locator("button.result-nav-btn", has_text="Combined table").first.click()
    page.wait_for_timeout(300)
    headers = page.locator("table.data-table:visible thead th").all_inner_texts()
    check(all(f"{o}: p-value" in headers for o in ["Genomics", "Transcriptomics", "Proteomics", "Metabolomics"]),
          "the table carries every omics' own p-value column")

    print("=== plain Fisher stays available and says what it ignores ===")
    page.locator("button", has_text="Start new run").first.click()
    page.wait_for_timeout(500)
    page.locator("select.text-input").filter(has=page.locator("option", has_text="Brown")).first.select_option("fisher")
    page.locator("button.ora-run-btn").click()
    for _ in range(60):
        if page.evaluate("STATE.aaPage") == "multiomics_results":
            break
        page.wait_for_timeout(500)
    ok = wait_for(page, "!!STATE.multiomics.run.results", 1500, "Fisher re-run")
    check(ok and page.evaluate("STATE.multiomics.run.results.blocks.KEGG.method") == "fisher", "the Fisher run completed")
    page.wait_for_timeout(400)
    page.locator("button.result-nav-btn", has_text="Overview").first.click()
    page.wait_for_timeout(300)
    body = page.inner_text("body")
    check("Shown for information only - Fisher does not use it" in body,
          "a Fisher run still shows the correlation table, marked as information only")
    check("as if independent" in body and "combined with Brown's method" not in body,
          "...and carries the independence caveat instead of Brown's")
    check(len(errors) == 0, f"no page error while switching method (got {errors})")

    print("=== back out, one level at a time ===")
    page.locator("button", has_text="Back to Multi-omics").first.click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "multiomics_home" and page.locator(".aa-genomics-grid button.module-btn").count() == 6,
          "Back from the results returns to the six-box screen")
    page.locator("button", has_text="Back to Annotation Analytics").first.click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "overview" and page.locator("button.module-btn.multiomics").count() == 1,
          "Back from the tile screen returns to the Annotation Analytics grid")

    print("=== console error check ===")
    check(len(errors) == 0, f"no console/page errors during the whole flow (got {errors})")
    browser.close()

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
