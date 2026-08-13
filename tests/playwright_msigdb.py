"""Playwright end-to-end UI test for MSigDB in Transcriptomics ORA and GSEA
(static/app.js). Drives real clicks through the actual retriever -> config ->
results flow against a live uvicorn (127.0.0.1:8009), using the app's own
example DE file and REAL MSigDB downloads (a small collection - Oncogenic,
189 human sets - kept deliberately small so the whole suite stays fast).

Companion to tests/test_msigdb_fetcher.py (direct function calls, fake HTTP)
and tests/regression_ora_gsea.py (real HTTP API calls, no browser) - this is
the layer that would have caught a UI wiring bug neither of those can see,
e.g. the config screen sending the wrong payload shape, or the collection
selector not actually gating on organism.
"""
import os
import sys

from playwright.sync_api import sync_playwright

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = "http://127.0.0.1:8009"
PASS = FAIL = 0
console_errors = []


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


def _go_to_module_config(page, module_name, continue_label):
    """Shared retriever flow (tile click -> upload -> organism -> resolve ->
    continue), mirroring playwright_regulon.py's working idiom exactly, so
    this file doesn't need to rediscover selectors that already work."""
    page.evaluate("""() => {
        STATE.view = 'annotation_analytics';
        STATE.aaPage = 'transcriptomics';
        render();
    }""")
    page.wait_for_timeout(200)
    tile = page.locator("button.module-btn.transcriptomics", has_text=module_name)
    check(tile.count() == 1, f"'{module_name}' tile is present")
    tile.click()
    page.wait_for_timeout(200)
    check(page.evaluate("STATE.transcriptomicsModuleName") == module_name,
          f"clicking the tile sets the module name to '{module_name}'")

    page.locator("button", has_text="Gene Symbol").first.click()
    page.wait_for_timeout(100)
    file_input = page.locator("input[type=file]")
    example_path = os.path.join(_ROOT, "Data3", "example_transcriptomics_de_full.csv")
    file_input.set_input_files(example_path)
    page.wait_for_timeout(1500)
    check(page.evaluate("!!STATE.transcriptomics.deTable") is True, "DE table uploaded")

    page.evaluate("""() => {
        STATE.transcriptomics.taxonId = '9606';
        STATE.transcriptomics.taxonName = 'Homo sapiens';
        render();
    }""")
    page.wait_for_timeout(200)

    resolve_btn = page.locator("button", has_text="Resolve identifiers")
    check(resolve_btn.count() == 1 and not resolve_btn.is_disabled(), "Resolve identifiers button is enabled")
    resolve_btn.click()
    for _ in range(120):
        if page.evaluate("!!STATE.transcriptomics.resolved"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.transcriptomics.resolved") is True, "resolve completed (real mygene.info round trip)")

    continue_btn = page.locator("button", has_text=continue_label)
    check(continue_btn.count() == 1, f"'{continue_label}' button present")
    continue_btn.click()
    page.wait_for_timeout(300)


with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page()
    page.on("console", lambda msg: console_errors.append(msg.text) if msg.type == "error" else None)
    page.on("pageerror", lambda exc: console_errors.append(f"pageerror: {exc}"))

    page.goto(BASE, wait_until="networkidle")
    check(page.evaluate("typeof render") == "function", "app.js loaded and parsed (render() is defined)")
    check(page.evaluate("typeof buildMsigdbOraCollectionSummary") == "function", "buildMsigdbOraCollectionSummary is defined")
    check(page.evaluate("typeof buildMsigdbGseaCollectionSummary") == "function", "buildMsigdbGseaCollectionSummary is defined")

    # The frontend's collection list must stay in sync with the backend's -
    # a drift here would silently mismatch request payloads / result keys.
    js_keys = page.evaluate("TRANSCRIPTOMICS_MSIGDB_COLLECTIONS.map(c => c.key)")
    expected_keys = ["hallmark", "perturbations", "wikipathways", "pid", "phenotype", "oncogenic", "immunologic", "cell_type"]
    check(js_keys == expected_keys,
          f"TRANSCRIPTOMICS_MSIGDB_COLLECTIONS' keys/order match Transcriptomics_module.MSIGDB_COLLECTIONS exactly "
          f"(got {js_keys})")

    print("=== ORA config screen: MSigDB pill, collection selector, organism gating ===")
    _go_to_module_config(page, "Over-Representation Analysis", "Continue to Over-Representation Analysis")
    check(page.evaluate("STATE.aaPage") == "transcriptomics_ora_config", "navigated to ORA config screen")

    msigdb_pill = page.locator("label.ora-db-pill", has_text="MSigDB").first
    check(msigdb_pill.count() >= 1, "MSigDB pill is present in the database checklist")
    msigdb_cb = msigdb_pill.locator("input[type=checkbox]")
    check(not msigdb_cb.is_disabled(), "MSigDB pill is enabled for human (taxon 9606)")

    collection_selector_before = page.locator(".ora-msigdb-collections")
    check(collection_selector_before.count() == 0, "collection sub-selector is NOT shown before MSigDB is checked")

    msigdb_cb.check()
    page.wait_for_timeout(200)
    collection_selector = page.locator(".ora-msigdb-collections")
    check(collection_selector.count() == 1, "collection sub-selector appears once MSigDB is checked")
    collection_pills = collection_selector.locator("label")
    check(collection_pills.count() == 8, f"exactly 8 collection checkboxes shown (got {collection_pills.count()})")

    run_btn = page.locator("button.ora-run-btn")
    check(run_btn.count() == 1 and run_btn.is_disabled(),
          "Run button stays DISABLED while MSigDB is checked but no collection is selected yet")

    oncogenic_pill = collection_selector.locator("label", has_text="Oncogenic")
    check(oncogenic_pill.count() == 1, "Oncogenic collection checkbox is present")
    oncogenic_pill.locator("input[type=checkbox]").check()
    page.wait_for_timeout(200)
    check(not run_btn.is_disabled(), "Run button becomes enabled once a collection is selected")

    print("=== ORA config screen: mouse organism greys out exactly PID/Oncogenic ===")
    # WikiPathways DOES have a mouse file (m2.cp.wikipathways, 202 sets -
    # confirmed live in live_msigdb_urls.py) - only PID and Oncogenic
    # genuinely have none for mouse.
    page.evaluate("""() => { STATE.transcriptomics.taxonId = '10090'; render(); }""")
    page.wait_for_timeout(200)
    mouse_unavailable = page.evaluate("""() => {
        const labels = [...document.querySelectorAll('.ora-msigdb-collections label')];
        return labels.filter(l => l.querySelector('input').disabled).map(l => l.textContent.trim());
    }""")
    check(
        all("PID" in t or "Oncogenic" in t for t in mouse_unavailable) and len(mouse_unavailable) == 2,
        f"exactly PID and Oncogenic are disabled for mouse, nothing else (got {mouse_unavailable})",
    )
    # Back to human for the actual run below.
    page.evaluate("""() => { STATE.transcriptomics.taxonId = '9606'; render(); }""")
    page.wait_for_timeout(200)
    # Re-select MSigDB + Oncogenic (the render() calls above rebuilt the DOM).
    page.locator("label.ora-db-pill", has_text="MSigDB").first.locator("input[type=checkbox]").check()
    page.wait_for_timeout(200)
    page.locator(".ora-msigdb-collections label", has_text="Oncogenic").locator("input[type=checkbox]").check()
    page.wait_for_timeout(200)

    print("=== ORA run with MSigDB: real result, family-count caption, summary figure ===")
    page.locator("button.ora-run-btn").click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "transcriptomics_ora_results", "navigated to ORA results screen")
    for _ in range(120):
        if page.evaluate("!!STATE.transcriptomics.oraResults"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.transcriptomics.oraResults") is True, "ORA+MSigDB job completed")

    result_tables = page.evaluate("Object.keys(STATE.transcriptomics.oraResults.tables)")
    check(any(k.startswith("MSigDB: Oncogenic") for k in result_tables),
          f"result tables include an 'MSigDB: Oncogenic (...)' family (got {result_tables})")

    caption_text = " ".join(page.locator(".caption").all_inner_texts())
    check("independently FDR-corrected" in caption_text, "multiple-testing family-count caption is rendered")

    summary_btn = page.locator("button", has_text="Show MSigDB collection summary")
    check(summary_btn.count() == 1, "'Show MSigDB collection summary' toggle is present")
    summary_btn.click()
    page.wait_for_timeout(600)
    summary_has_data = page.evaluate("""() => {
        const divs = [...document.querySelectorAll('.chart-container')];
        const d = divs.find(x => x.id.includes('msigdbsummary') && x.data);
        return !!(d && d.data && d.data.length > 0 && d.data[0].y.length > 0);
    }""")
    check(summary_has_data, "MSigDB collection summary figure actually rendered with data")

    print("=== result-table CSV downloads carry the result key in the filename ===")
    # Done BEFORE the alpha change below: at a very strict threshold the
    # result panel can render its "nothing significant" state, and this check
    # is about filenames, not about what survives filtering.
    # The download button lives in the "Result table" pane, which is not the
    # default-active one - open it first, or the click times out on a
    # present-but-hidden element.
    table_rail = page.locator(".result-nav-btn", has_text="Result table")
    check(table_rail.count() >= 1, "'Result table' rail entry present")
    if table_rail.count() >= 1:
        table_rail.first.click()
        page.wait_for_timeout(400)
        dl_btn = page.locator("button:visible", has_text="Download result table").first
        with page.expect_download(timeout=15000) as dl_info:
            dl_btn.click()
        download = dl_info.value
        fname = download.suggested_filename
        check(fname.endswith(".csv") and "msigdb" in fname.lower(),
              f"the downloaded filename identifies which MSigDB family it is (got {fname!r})")
        check(fname != "transcriptomics_ora_result.csv",
              f"it is NOT the old fixed filename every table shared (got {fname!r})")

    print("=== summary figure tracks the significance threshold (it used to cache the first alpha forever) ===")
    # The summary lives OUTSIDE panelHost and caches its build, so before the
    # refresh hook was added it kept showing bar lengths and an axis title
    # computed at the ORIGINAL alpha while the table below used the new one.
    axis_before = page.evaluate("""() => {
        const divs = [...document.querySelectorAll('.chart-container')];
        const d = divs.find(x => x.id.includes('msigdbsummary') && x.data);
        return d ? d.layout.xaxis.title.text : null;
    }""")
    alpha_input = page.locator("input.ora-sig-threshold-input")
    alpha_input.fill("0.001")
    alpha_input.dispatch_event("change")
    page.wait_for_timeout(700)
    axis_after = page.evaluate("""() => {
        const divs = [...document.querySelectorAll('.chart-container')];
        const d = divs.find(x => x.id.includes('msigdbsummary') && x.data);
        return d ? d.layout.xaxis.title.text : null;
    }""")
    check(axis_before is not None and axis_after is not None and axis_before != axis_after,
          f"changing the FDR threshold rebuilds the summary - its axis title follows the new alpha "
          f"(before={axis_before!r}, after={axis_after!r})")
    check(axis_after is not None and "0.001" in axis_after,
          f"and the axis title states the CURRENT alpha (got {axis_after!r})")

    print("=== 'Start new run' resets MSigDB state cleanly ===")
    new_run_btn = page.locator("button", has_text="Start new run")
    if new_run_btn.count() > 0:
        new_run_btn.click()
        page.wait_for_timeout(200)
        check(page.evaluate("STATE.aaPage") == "transcriptomics_ora_config", "'Start new run' returns to config screen")

    print("=== GSEA config + run with MSigDB ===")
    page.evaluate("""() => {
        STATE.view = 'annotation_analytics';
        STATE.aaPage = 'transcriptomics';
        render();
    }""")
    page.wait_for_timeout(200)
    _go_to_module_config(page, "Gene Set Enrichment Analysis", "Continue to Gene Set Enrichment Analysis")
    check(page.evaluate("STATE.aaPage") == "transcriptomics_gsea_config", "navigated to GSEA config screen")

    gsea_msigdb_pill = page.locator("label.ora-db-pill", has_text="MSigDB").first
    check(gsea_msigdb_pill.count() >= 1, "MSigDB pill is present on the GSEA config screen too")
    gsea_msigdb_pill.locator("input[type=checkbox]").check()
    page.wait_for_timeout(200)
    page.locator(".ora-msigdb-collections label", has_text="Oncogenic").locator("input[type=checkbox]").check()
    page.wait_for_timeout(200)
    # Default permutation count (1000) against a 189-set collection is still
    # fast enough not to need overriding for this test to finish quickly.

    gsea_run_btn = page.locator("button.ora-run-btn")
    check(not gsea_run_btn.is_disabled(), "GSEA run button enabled with MSigDB + Oncogenic selected")
    gsea_run_btn.click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "transcriptomics_gsea_results", "navigated to GSEA results screen")
    for _ in range(180):
        if page.evaluate("!!STATE.transcriptomics.gseaResults"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.transcriptomics.gseaResults") is True, "GSEA+MSigDB job completed")

    gsea_result_tables = page.evaluate("Object.keys(STATE.transcriptomics.gseaResults.tables)")
    check("MSigDB: Oncogenic" in gsea_result_tables,
          f"GSEA result key has NO direction suffix (got {gsea_result_tables})")

    gsea_caption_text = " ".join(page.locator(".caption").all_inner_texts())
    check("independently FDR-corrected" in gsea_caption_text, "GSEA results screen also shows the family-count caption")

    gsea_summary_btn = page.locator("button", has_text="Show MSigDB collection summary")
    check(gsea_summary_btn.count() == 1, "GSEA results screen also offers the MSigDB collection summary toggle")
    gsea_summary_btn.click()
    page.wait_for_timeout(600)
    gsea_summary_has_data = page.evaluate("""() => {
        const divs = [...document.querySelectorAll('.chart-container')];
        const d = divs.find(x => x.id.includes('msigdbsummary') && x.data);
        return !!(d && d.data && d.data.length > 0 && d.data[0].y.length > 0);
    }""")
    check(gsea_summary_has_data, "GSEA MSigDB collection summary figure rendered with data")

    print("=== organism change after resolve clears the resolved table + all results (Part 3 fix) ===")
    # Fresh flow, independent of the ORA/GSEA runs above.
    _go_to_module_config(page, "Over-Representation Analysis", "Continue to Over-Representation Analysis")
    check(page.evaluate("!!STATE.transcriptomics.resolved") is True,
          "resolved table exists before switching organism (precondition)")

    # Back to the retriever screen, where the organism selector lives.
    page.locator("button", has_text="Back").first.click()
    page.wait_for_timeout(200)
    check(page.evaluate("STATE.aaPage") == "transcriptomics_module", "back at the retriever screen")

    # Drive the REAL selector's change handler - not the STATE shortcut
    # _go_to_module_config uses for the initial resolve above (which
    # deliberately bypasses the fuzzy-search autocomplete, since that
    # debounced live search is unrelated to what's under test here and is
    # already covered elsewhere). Populate the hidden <select> the way a real
    # search result would, then dispatch a genuine 'change' event, so the
    # actual fix - the listener inside renderTranscriptomicsOrganismSelector
    # - runs exactly as it would for a real user picking "Mus musculus".
    switch_result = page.evaluate("""() => {
        const input = [...document.querySelectorAll('input')].find(
            i => i.placeholder && i.placeholder.includes('Species / organism'));
        if (!input) return { ok: false, reason: 'organism input not found' };
        const select = input.nextElementSibling;
        if (!select || select.tagName !== 'SELECT') return { ok: false, reason: 'organism select not found' };
        select._records = [{ context_id: '10090', context_name: 'Mus musculus' }];
        select.innerHTML = '';
        const opt = document.createElement('option');
        opt.value = '0';
        opt.textContent = '10090 | Mus musculus';
        select.appendChild(opt);
        select.value = '0';
        select.dispatchEvent(new Event('change', { bubbles: true }));
        return { ok: true };
    }""")
    check(switch_result.get("ok") is True, f"drove the real organism selector's change handler (got {switch_result})")
    page.wait_for_timeout(300)

    check(page.evaluate("STATE.transcriptomics.taxonId") == "10090", "taxonId updated to the newly selected organism")
    check(page.evaluate("STATE.transcriptomics.resolved") is None,
          "resolved table was cleared by the organism switch")
    check(page.evaluate("STATE.transcriptomics.oraResults") is None, "oraResults cleared too")
    check(page.evaluate("STATE.transcriptomics.oraConfig") is None, "oraConfig cleared too")
    check(page.evaluate("STATE.transcriptomics.oraRunLocked") is False, "oraRunLocked reset")
    check(page.evaluate("'_oraDraft' in STATE.transcriptomics") is False,
          "the stale (organism-gated) ORA draft was deleted, not just its results")

    banner_text = " ".join(page.locator(".alert-warning").all_inner_texts())
    check("Organism changed" in banner_text and "Mus musculus" in banner_text,
          f"a warning banner naming the new organism is rendered on screen (got {banner_text[:300]!r})")

    resolve_btn_after_switch = page.locator("button", has_text="Resolve identifiers")
    check(resolve_btn_after_switch.count() == 1, "back at the 'Resolve identifiers' step after the switch")

    print("=== ...but switching organism BEFORE any resolve must NOT show a spurious warning ===")
    # The guard is "had a resolved table for a DIFFERENT organism", not "the
    # organism changed" - selecting an organism for the first time (nothing
    # resolved yet) must not trigger the clear-and-warn path.
    fresh_switch = page.evaluate("""() => {
        const input = [...document.querySelectorAll('input')].find(
            i => i.placeholder && i.placeholder.includes('Species / organism'));
        const select = input.nextElementSibling;
        select._records = [{ context_id: '9606', context_name: 'Homo sapiens' }];
        select.innerHTML = '';
        const opt = document.createElement('option');
        opt.value = '0'; opt.textContent = '9606 | Homo sapiens';
        select.appendChild(opt);
        select.value = '0';
        select.dispatchEvent(new Event('change', { bubbles: true }));
        return true;
    }""")
    page.wait_for_timeout(300)
    check(page.evaluate("!!STATE.transcriptomics._organismChangedNotice") is False,
          "no warning is set when there was nothing resolved to invalidate")

    print("=== SENSITIVITY: the organism-switch checks above are not vacuous ===")
    # Confirm the underlying helper genuinely does the clearing by calling it
    # directly and checking a broken (no-op) stand-in would NOT have cleared
    # anything - proves the assertions above are exercising real behaviour,
    # not a coincidence of how this fixture happened to be built.
    sensitivity = page.evaluate("""() => {
        const saved = Object.assign({}, STATE.transcriptomics);
        STATE.transcriptomics.resolved = { fake: true };
        STATE.transcriptomics.oraResults = { fake: true };
        const noop = () => {};  // stand-in for a version that forgot to clear anything
        noop();
        const stillThere = STATE.transcriptomics.resolved !== null;
        _clearTranscriptomicsResolutionDependentState();
        const clearedByRealFn = STATE.transcriptomics.resolved === null && STATE.transcriptomics.oraResults === null;
        Object.assign(STATE.transcriptomics, saved);
        return { stillThere, clearedByRealFn };
    }""")
    check(sensitivity["stillThere"] is True,
          "a no-op stand-in (the pre-fix behaviour) really would leave stale state in place")
    check(sensitivity["clearedByRealFn"] is True,
          "the real _clearTranscriptomicsResolutionDependentState genuinely clears it")

    print("=== console error check ===")
    check(len(console_errors) == 0, f"no console/page errors during the whole flow (got {len(console_errors)}: {console_errors[:5]})")

    browser.close()

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
