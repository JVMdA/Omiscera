"""FIRST real-browser coverage of the Proteomics module (this whole
submodule is new - Milestone 1: foundation + Over-Representation Analysis).

Drives the REAL flow end to end: tile -> organism selector -> upload the
bundled example -> resolve (real UniProt ID Mapping) -> configure -> run
(real GO/InterPro/Complex Portal/UniProt keyword fetches) -> results, then
separately re-uploads the deliberately messy MaxQuant-style example to
verify contaminant/protein-group handling is visible on screen, not just
correct server-side (already unit-tested in test_proteomics_retriever.py).

Needs a live server and real network access - not mocked, takes a few
minutes (GO's QuickGO fetch + hierarchy propagation is the slow part).
"""
import os

from playwright.sync_api import sync_playwright

BASE = "http://127.0.0.1:8009"
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
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

    print("=== the Proteomics tile is present and enabled ===")
    page.evaluate("""() => { STATE.view = 'annotation_analytics'; STATE.aaPage = 'overview'; render(); }""")
    page.wait_for_timeout(200)
    proteomics_tile = page.locator("button.module-btn.proteomics")
    check(proteomics_tile.count() >= 1, "top-level 'Proteomics' tile is present")
    check(not proteomics_tile.first.is_disabled(), "the tile is ENABLED (no longer a dead button)")
    proteomics_tile.first.click()
    page.wait_for_timeout(200)
    check(page.evaluate("STATE.aaPage") == "proteomics", "clicking the tile navigated to the Proteomics grid")

    print("=== ORA sub-tile enabled, other three disabled ===")
    ora_btn = page.locator("button.module-btn.proteomics", has_text="Over-Representation Analysis")
    check(ora_btn.count() >= 1 and not ora_btn.first.is_disabled(), "ORA sub-tile present and enabled")
    for name in ["Quantitative Enrichment", "PPI Network / Module Detection", "PTM / Phosphosite Enrichment"]:
        btn = page.locator("button.module-btn.proteomics", has_text=name)
        check(btn.count() >= 1 and btn.first.is_disabled(), f"'{name}' correctly disabled (coming soon)")
    ora_btn.first.click()
    page.wait_for_timeout(200)
    check(page.evaluate("STATE.aaPage") == "proteomics_module", "navigated to the upload/resolve screen")

    print("=== organism + upload + resolve (real UniProt ID Mapping) ===")
    page.evaluate("""() => { STATE.proteomics.taxonId = '9606'; STATE.proteomics.taxonName = 'Homo sapiens'; render(); }""")
    page.wait_for_timeout(150)
    page.locator("input[type=file]").set_input_files(
        os.path.join(_ROOT, "Data3", "example_proteomics_de_full.csv")
    )
    page.wait_for_timeout(1500)
    check(page.evaluate("!!STATE.proteomics.deTable"), "upload succeeded, preview table populated")

    resolve_btn = page.locator("button", has_text="Resolve accessions")
    check(resolve_btn.count() >= 1 and not resolve_btn.first.is_disabled(), "Resolve accessions button enabled")
    resolve_btn.first.click()
    for _ in range(180):
        if page.evaluate("!!STATE.proteomics.resolved"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.proteomics.resolved") is True, "resolve completed (real UniProt ID Mapping round trip)")
    meta = page.evaluate("STATE.proteomics.resolved.metrics")
    check(meta["n_mapped"] > 0, f"real accessions resolved (n_mapped={meta['n_mapped']})")
    check(meta["n_unmapped"] == 0, f"the clean example has zero unmapped accessions (got {meta['n_unmapped']})")

    print("=== continue -> configure -> run (fast databases only, for speed) ===")
    page.locator("button", has_text="Continue to Over-Representation Analysis").click()
    page.wait_for_timeout(200)
    check(page.evaluate("STATE.aaPage") == "proteomics_ora_config", "navigated to the ORA config screen")

    for db_name in ["InterPro", "Complex Portal", "UniProt (keywords)"]:
        page.locator("label.ora-db-pill", has_text=db_name).locator("input[type=checkbox]").check()
    page.wait_for_timeout(100)
    run_btn = page.locator("button.ora-run-btn")
    check(not run_btn.is_disabled(), "Run button enabled with 3 databases selected")
    run_btn.click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "proteomics_ora_results", "navigated to the results screen")

    for _ in range(180):
        if page.evaluate("!!STATE.proteomics.oraResults"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.proteomics.oraResults") is True, "the ORA job completed")

    print("=== results: real, non-empty, biologically sensible hits ===")
    summary = page.evaluate("STATE.proteomics.oraResults.summary")
    check(summary["n_sig_all"] > 0, f"the run found significant proteins (n_sig_all={summary['n_sig_all']})")
    tables = page.evaluate("STATE.proteomics.oraResults.tables")
    cpx_all = tables.get("Complex Portal (all)")
    check(cpx_all is not None and cpx_all["n_terms_tested"] > 0,
          "Complex Portal (all) table exists and has tested terms")
    if cpx_all and cpx_all["table"]["rows"]:
        top_name = cpx_all["table"]["rows"][0].get("Term Name", "")
        check("ribosom" in top_name.lower() or "proteasome" in top_name.lower() or "spliceosom" in top_name.lower(),
              f"Complex Portal's top hit is one of the real complexes seeded into the example (got: {top_name!r})")

    check(page.locator("select.transcriptomics-ora-result-select").count() >= 1, "the result-table selector is present")
    check(page.locator(".ora-db-card").count() >= 1, "a results card rendered")

    print("=== figures render on the current panel ===")
    # "Show dot plot" lives inside the "Protein set overview" result-nav pane,
    # which is not the default-active one (that's "Overview") - select it
    # first, or the button exists in the DOM but stays display:none.
    overview_rail_btn = page.locator("button.result-nav-btn", has_text="Protein set overview")
    check(overview_rail_btn.count() >= 1, "'Protein set overview' rail entry is present (this table has significant rows)")
    if overview_rail_btn.count() >= 1:
        overview_rail_btn.first.click()
        page.wait_for_timeout(200)
        dot_btn = page.locator("button.chart-toggle-btn", has_text="Show dot plot")
        check(dot_btn.count() >= 1, "'Show dot plot' button present after selecting the pane")
        dot_btn.first.click()
        page.wait_for_timeout(600)
        check(page.locator(".js-plotly-plot").count() >= 1, "dot plot rendered a real Plotly figure")

    print("=== the messy MaxQuant example: contaminant filtering visible on screen ===")
    page.evaluate("""() => { STATE.aaPage = 'proteomics'; STATE.proteomicsModuleName = ''; render(); }""")
    page.wait_for_timeout(150)
    page.locator("button.module-btn.proteomics", has_text="Over-Representation Analysis").click()
    page.wait_for_timeout(150)
    page.evaluate("""() => { STATE.proteomics.taxonId = '9606'; STATE.proteomics.taxonName = 'Homo sapiens'; render(); }""")
    page.wait_for_timeout(150)
    page.locator("input[type=file]").set_input_files(
        os.path.join(_ROOT, "Data3", "example_proteomics_maxquant_messy.csv")
    )
    page.wait_for_timeout(1000)
    page.locator("button", has_text="Resolve accessions").first.click()
    for _ in range(120):
        if page.evaluate("!!STATE.proteomics.resolved"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.proteomics.resolved") is True, "messy example: resolve completed")
    meta2 = page.evaluate("STATE.proteomics.resolved.metrics")
    n_removed = meta2["contaminants"]["n_total_removed"]
    check(n_removed >= 3, f"messy example: contaminant/decoy rows were actually removed (n_total_removed={n_removed})")
    page.wait_for_timeout(200)
    caveat_visible = page.evaluate("""() => {
        const boxes = Array.from(document.querySelectorAll('.alert-box, .alert'));
        return boxes.some(b => b.textContent.includes('contaminant'));
    }""")
    check(caveat_visible, "the contaminant-removal count is shown on screen, not just in the payload")

    print("=== console error check ===")
    check(len(errors) == 0, f"no console/page errors during the whole flow (got {errors})")

    browser.close()

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
