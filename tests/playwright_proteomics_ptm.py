"""Real-browser coverage of the Proteomics "PTM / Phosphosite Enrichment"
tile, the 4th and final Proteomics tile - newly enabled (previously in
PROTEOMICS_TILES_NOT_YET_IMPLEMENTED).

Drives the REAL flow end to end: tile enabled -> UniProt-taxonomy organism
selector -> upload the bundled site-level example -> resolve (real UniProt
ID Mapping, isoform-preserving) -> continue to config -> run (real iPTMnet
fetch + kinase activity ULM + kinase substrate ORA) -> results, checking that
MAPK1 (the kinase whose substrate sites were simulated strongly UP in the
bundled example - see Data3/example_proteomics_phosphosites.csv) is surfaced
correctly in both the activity table and the substrate-enrichment table, and
that the CSV download works.

Needs a live server and real network access - not mocked. The organism
search hits a live UniProt taxonomy lookup (can take longer than STRING's
local species file), hence the generous wait after typing the query.
"""
import os

from playwright.sync_api import sync_playwright

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


with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page()
    errors = []
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    page.goto(BASE, wait_until="networkidle")
    page.wait_for_timeout(300)
    check(len(errors) == 0, f"app.js parsed with no top-level error (got {errors})")

    page.evaluate("""() => { STATE.view = 'annotation_analytics'; STATE.aaPage = 'overview'; render(); }""")
    page.wait_for_timeout(150)
    page.locator("button.module-btn.proteomics").first.click()
    page.wait_for_timeout(150)

    print("=== PTM tile is enabled (no longer 'coming soon') ===")
    ptm_btn = page.locator("button.module-btn.proteomics", has_text="PTM / Phosphosite Enrichment")
    check(ptm_btn.count() >= 1, "PTM tile button present")
    check(not ptm_btn.first.is_disabled(), "PTM tile is ENABLED")

    ptm_btn.first.click()
    page.wait_for_timeout(150)
    check(page.evaluate("STATE.aaPage") == "proteomics_module", "clicking the tile navigated to proteomics_module")
    check(page.evaluate("STATE.proteomicsModuleName") == "PTM / Phosphosite Enrichment", "proteomicsModuleName set correctly")

    body_text = page.inner_text("body")
    check("Upload phosphosite differential-abundance table" in body_text,
          "site-level upload section is shown (not the protein-level one)")

    print("=== organism selection (real UniProt taxonomy search) ===")
    org_input = page.locator("input.text-input").first
    org_input.fill("Homo sapiens")
    page.wait_for_timeout(2500)
    options = page.locator("select").first.locator("option").all_inner_texts()
    idx = next((i for i, o in enumerate(options) if "9606" in o), None)
    check(idx is not None, f"organism search for 'Homo sapiens' returns taxon 9606 (got {options[:5]})")
    page.locator("select").first.select_option(index=idx)
    page.wait_for_timeout(150)
    check(page.evaluate("STATE.proteomics.taxonId") == "9606", "organism selection sets taxonId")

    print("=== upload the bundled PTM example, resolve, reach the config screen ===")
    import requests
    ex_bytes = requests.get(BASE + "/api/proteomics/ptm/example").content
    tmp_csv = os.path.join(os.environ.get("TEMP", "."), "pw_ptm_example.csv")
    with open(tmp_csv, "wb") as fh:
        fh.write(ex_bytes)
    with page.expect_file_chooser() as fc_info:
        page.locator("text=Choose CSV file").click()
    fc_info.value.set_files(tmp_csv)
    page.wait_for_timeout(1500)
    check(page.evaluate("!!STATE.proteomics.ptmSitesTable"), "example file uploaded and previewed")

    resolve_btn = page.locator("button", has_text="Resolve sites")
    check(resolve_btn.count() >= 1 and not resolve_btn.first.is_disabled(), "Resolve sites button enabled")
    resolve_btn.first.click()
    page.wait_for_function("!!STATE.proteomics.resolvedSites", timeout=90000)
    check(True, "resolve completed")

    continue_btn = page.locator("button", has_text="Continue to PTM / Phosphosite Enrichment")
    check(continue_btn.count() >= 1, "Continue button present after resolve")
    continue_btn.first.click()
    page.wait_for_timeout(200)
    check(page.evaluate("STATE.aaPage") == "proteomics_ptm_config", "navigated to proteomics_ptm_config")

    body_text = page.inner_text("body")
    check("Proteomics — PTM / Phosphosite Enrichment: setup" in body_text, "config screen title correct")

    print("=== run the analysis for real (500 permutations) ===")
    perm_selects = page.locator("select.chart-option-select")
    picked = False
    for i in range(perm_selects.count()):
        opts = perm_selects.nth(i).locator("option").all_inner_texts()
        if opts == ["500", "1000", "2000", "5000"]:
            perm_selects.nth(i).select_option("500")
            picked = True
            break
    check(picked, "found and set the permutation-count selector to 500")

    run_btn = page.locator("button", has_text="Run").first
    check(not run_btn.is_disabled(), "Run button enabled")
    run_btn.click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "proteomics_ptm_results", "navigated to results screen")

    page.wait_for_function("!!STATE.proteomics.ptmActivityResults", timeout=180000)
    check(True, "run completed (activity results populated)")

    print("=== results content ===")
    body_text = page.inner_text("body")
    check("Kinases tested" in body_text, "overview shows 'Kinases tested'")

    activity_tab = page.locator("button.result-nav-btn", has_text="Kinase activity")
    check(activity_tab.count() >= 1, "'Kinase activity' nav tab exists")
    activity_tab.first.click()
    page.wait_for_timeout(200)
    activity_text = page.inner_text("body")
    check("MAPK1" in activity_text, "MAPK1 (the simulated-active kinase) appears in the activity results")
    check(page.locator("label.chart-option-control", has_text="Kinases shown (bar plot)").count() >= 1,
          "the kinase activity bar plot has a 'Kinases shown' topN control (bug found 2026-09-21: "
          "this used to be hardcoded with no control at all)")

    ora_tab = page.locator("button.result-nav-btn", has_text="Kinase substrate enrichment")
    check(ora_tab.count() >= 1, "'Kinase substrate enrichment' nav tab exists")
    ora_tab.first.click()
    page.wait_for_timeout(200)
    ora_text = page.inner_text("body")
    check("MAPK1" in ora_text, "MAPK1 appears in the substrate enrichment results too")

    print("=== dot/bar plot chart controls (bug found 2026-09-21: this pane used to hardcode "
          "topN=20/fontSize=12/default colors with NO controls at all) ===")
    for label in ["Kinases shown", "Rank & show by", "Font size", "Dot gradient", "Bar color"]:
        cnt = page.locator("label.chart-option-control", has_text=label).count()
        check(cnt >= 1, f"'{label}' control is present on the Kinase substrate enrichment pane (got {cnt})")

    dot_btn = page.locator("button", has_text="Show dot plot")
    dot_btn.first.click()
    page.wait_for_timeout(1000)
    check(page.locator(".chart-container:visible .plotly").count() >= 1,
          "the dot plot actually renders (a visible Plotly chart exists)")

    gradient_sel = page.locator("label.chart-option-control", has_text="Dot gradient").locator("select")
    grad_options = gradient_sel.locator("option").all_inner_texts()
    check(len(grad_options) > 1, f"Dot gradient offers multiple choices (got {grad_options})")
    other_grad = next((o for o in grad_options if o != gradient_sel.input_value()), grad_options[-1])
    gradient_sel.select_option(label=other_grad)
    page.wait_for_timeout(500)
    check(len(errors) == 0, f"changing the dot gradient redraws without a console error (got {errors})")

    fontsize_sel = page.locator("label.chart-option-control", has_text="Font size").locator("select")
    fontsize_sel.select_option("18")
    page.wait_for_timeout(500)
    check(len(errors) == 0, f"changing font size redraws without a console error (got {errors})")

    topn_sel = page.locator(
        'label.chart-option-control:has(span.chart-option-label:text-is("Kinases shown"))').locator("select")
    topn_sel.select_option("5")
    page.wait_for_timeout(500)
    check(len(errors) == 0, f"changing 'Kinases shown' redraws without a console error (got {errors})")

    bar_btn = page.locator('button:text-is("Show bar plot")')
    bar_btn.first.click()
    page.wait_for_timeout(1000)
    barcolor_sel = page.locator("label.chart-option-control", has_text="Bar color").locator("select")
    bar_options = barcolor_sel.locator("option").all_inner_texts()
    other_bar = next((o for o in bar_options if o != barcolor_sel.input_value()), bar_options[-1])
    barcolor_sel.select_option(label=other_bar)
    page.wait_for_timeout(500)
    check(len(errors) == 0, f"changing bar color redraws without a console error (got {errors})")

    print("=== CSV download ===")
    activity_tab.first.click()
    page.wait_for_timeout(150)
    dl_btn = page.locator("button", has_text="Download kinase activity (CSV)")
    check(dl_btn.count() >= 1, "activity CSV download button present")
    with page.expect_download() as dl_info:
        dl_btn.first.click()
    dl = dl_info.value
    check("proteomics_ptm" in dl.suggested_filename, f"downloaded filename is PTM-specific (got {dl.suggested_filename})")

    print("=== console error check ===")
    check(len(errors) == 0, f"no console/page errors during the whole flow (got {errors})")

    browser.close()

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
