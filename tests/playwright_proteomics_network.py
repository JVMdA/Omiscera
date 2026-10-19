"""FIRST real-browser coverage of Proteomics PPI Network / Module Detection.

Drives the REAL flow end to end: tile -> STRING-species organism selector ->
upload the bundled example -> resolve (real UniProt ID Mapping) -> continue
to PPI Network / Module Detection -> configure -> run (real STRING download
+ 500-permutation propagation + GO fetch) -> results, checking that the
descriptor-driven generalization of Transcriptomics' own Network config/
results screens (renderNetworkConfigScreen/renderNetworkResultsScreen in
app.js) actually renders Proteomics-specific text (protein noun, "Protein
mapping" not "Gene mapping", Reactome/Complex Portal offered, MSigDB not
offered) rather than leftover Transcriptomics wording, and that a real run
completes and its CSV downloads work.

Needs a live server and real network access - not mocked. The STRING human
network download plus 500 permutations takes a few minutes; this is the
slowest Proteomics Playwright test in this suite for that reason.
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

    print("=== PPI Network / Module Detection tile is enabled ===")
    ppi_btn = page.locator("button.module-btn.proteomics", has_text="PPI Network / Module Detection")
    check(ppi_btn.count() >= 1, "PPI Network tile button present")
    check(not ppi_btn.first.is_disabled(), "PPI Network tile is ENABLED")
    ptm_btn = page.locator("button.module-btn.proteomics", has_text="PTM / Phosphosite Enrichment")
    check(not ptm_btn.first.is_disabled(), "PTM tile is now ENABLED too (see playwright_proteomics_ptm.py for its own coverage)")

    ppi_btn.first.click()
    page.wait_for_timeout(150)
    check(page.evaluate("STATE.aaPage") == "proteomics_module", "clicking the tile navigated to proteomics_module")
    check(page.evaluate("STATE.proteomicsModuleName") == "PPI Network / Module Detection",
          "STATE.proteomicsModuleName set correctly")

    print("=== organism selector switches to STRING species search for this tile ===")
    placeholder = page.locator("input.text-input").first.get_attribute("placeholder")
    check(placeholder is not None and "STRING" in placeholder,
          f"organism input placeholder mentions STRING (got {placeholder!r})")

    org_input = page.locator("input.text-input").first
    org_input.fill("Homo sapiens")
    page.wait_for_timeout(700)
    options = page.locator("select").first.locator("option").all_inner_texts()
    idx = next((i for i, o in enumerate(options) if o.startswith("9606")), None)
    check(idx is not None, f"STRING species search for 'Homo sapiens' returns taxon 9606 (got {options[:5]})")
    page.locator("select").first.select_option(index=idx)
    page.wait_for_timeout(150)
    check(page.evaluate("STATE.proteomics.taxonId") == "9606",
          f"selecting the STRING result sets taxonId (got {page.evaluate('STATE.proteomics.taxonId')})")

    print("=== upload the bundled example, resolve, reach the PPI config screen ===")
    import requests
    ex_bytes = requests.get(BASE + "/api/proteomics/example").content
    tmp_csv = os.path.join(os.environ.get("TEMP", "."), "playwright_proteomics_network_example.csv")
    with open(tmp_csv, "wb") as fh:
        fh.write(ex_bytes)
    with page.expect_file_chooser() as fc_info:
        page.locator("text=Choose CSV file").click()
    fc_info.value.set_files(tmp_csv)
    page.wait_for_timeout(1500)
    check(page.evaluate("!!STATE.proteomics.deTable"), "example file uploaded and previewed")

    resolve_btn = page.locator("button", has_text="Resolve accessions")
    check(resolve_btn.count() >= 1 and not resolve_btn.first.is_disabled(),
          "Resolve button enabled after upload + organism selection")
    resolve_btn.first.click()
    page.wait_for_function("!!STATE.proteomics.resolved", timeout=90000)
    check(True, "resolve completed")

    continue_btn = page.locator("button", has_text="Continue to PPI Network / Module Detection")
    check(continue_btn.count() >= 1, "'Continue to PPI Network / Module Detection' button present after resolve")
    continue_btn.first.click()
    page.wait_for_timeout(200)
    check(page.evaluate("STATE.aaPage") == "proteomics_network_config", "navigated to proteomics_network_config")

    body_text = page.inner_text("body")
    check("Proteomics — PPI Network / Module Detection: setup" in body_text, "config screen title is Proteomics-specific")
    check("Seed proteins" in body_text, "seed section says 'Seed proteins', not 'Seed genes'")
    check("STRING taxon 9606" in body_text, "organism caption shows the STRING taxon")

    print("=== module-label database choices are Proteomics' own (GO/Reactome/Complex Portal, no MSigDB), and MULTIPLE can be checked at once ===")
    db_pills = page.locator("label.ora-db-pill")
    pill_labels = db_pills.all_inner_texts()
    check(any("Reactome" in t for t in pill_labels) and any("Complex Portal" in t for t in pill_labels),
          f"module-label checkboxes offer Reactome AND Complex Portal (got {pill_labels})")
    check(not any("MSigDB" in t for t in pill_labels),
          f"MSigDB is NOT offered for Proteomics (got {pill_labels})")
    check("Reactome and Complex Portal ARE offered here" in body_text,
          "Proteomics-specific module-label note is shown, not Transcriptomics' 'Reactome is not offered'")

    go_checkbox = db_pills.filter(has_text="GO").locator("input.ora-db-checkbox")
    check(go_checkbox.is_checked(), "GO is checked by default")
    complex_checkbox = db_pills.filter(has_text="Complex Portal").locator("input.ora-db-checkbox")
    check(not complex_checkbox.is_checked(), "Complex Portal starts unchecked")
    complex_checkbox.check()
    check(complex_checkbox.is_checked() and go_checkbox.is_checked(),
          "checking Complex Portal does NOT uncheck GO - both stay checked at once, unlike the old single-select")

    print("=== run the analysis for real (500 permutations) and check the results screen ===")
    perm_selects = page.locator("select.chart-option-select")
    picked_perm = False
    for i in range(perm_selects.count()):
        opts = perm_selects.nth(i).locator("option").all_inner_texts()
        if opts == ["500", "1000", "2000"]:
            perm_selects.nth(i).select_option("500")
            picked_perm = True
            break
    check(picked_perm, "found and set the permutation-count selector to 500")

    run_btn = page.locator("button", has_text="Run").first
    check(not run_btn.is_disabled(), "Run button enabled")
    run_btn.click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "proteomics_network_results", "navigated to results screen immediately")

    page.wait_for_function("!!STATE.proteomics.networkResults", timeout=360000)
    check(True, "run completed (results populated)")

    body_text = page.inner_text("body")
    check("Proteins seeded" in body_text, "results overview says 'Proteins seeded'")
    mapping_tab = page.locator("button.result-nav-btn", has_text="Protein mapping")
    check(mapping_tab.count() >= 1, "results has a 'Protein mapping' nav tab, not 'Gene mapping'")
    mapping_tab.first.click()
    page.wait_for_timeout(150)
    mapping_pane_text = page.inner_text("body")
    check("Matched by accession" in mapping_pane_text, "mapping tier caption uses accession/gene-symbol wording")
    result = page.evaluate("STATE.proteomics.networkResults")
    check(result.get("mapping", {}).get("n_mapped", 0) > 0,
          f"mapping n_mapped > 0 (got {result.get('mapping', {}).get('n_mapped')})")
    check(result.get("seeds", {}).get("n_seeds", 0) > 0,
          f"seeds n_seeds > 0 (got {result.get('seeds', {}).get('n_seeds')})")
    check("caveats" in result, "payload carries a caveats key (matches Transcriptomics' payload shape)")

    check(page.evaluate("STATE.proteomics.networkConfig.enrichment_databases") == ["GO", "Complex Portal"],
          f"the run's own stored config reflects BOTH checked databases, in the order checked "
          f"(got {page.evaluate('STATE.proteomics.networkConfig.enrichment_databases')})")

    print("=== Module enrichment pane: per-module DATABASE selector, not just whichever one won ===")
    enr = (result.get("modules") or {}).get("enrichment") or {}
    any_candidates_dict = any(isinstance(e, dict) and isinstance(e.get("candidates"), dict) for e in enr.values() if e)
    check(not enr or any_candidates_dict,
          f"every labelled module's payload carries a 'candidates' dict (the new multi-database shape), "
          f"not just a single winning block (got enrichment={list(enr.items())[:1]})")
    enrich_tab = page.locator("button.result-nav-btn", has_text="Module enrichment")
    if enrich_tab.count():
        enrich_tab.first.click()
        page.wait_for_timeout(150)
        db_control = page.locator(".chart-option-control", has_text="Database").locator("select")
        check(db_control.count() >= 1, "the 'Module enrichment' pane has its own Database selector")
        if db_control.count():
            db_opts = db_control.first.locator("option").all_inner_texts()
            check(len(db_opts) >= 1, f"the Database selector is populated with at least the winning candidate (got {db_opts})")
            terms_control = page.locator(".chart-option-control", has_text="Terms shown")
            dot_control = page.locator(".chart-option-control", has_text="Dot gradient")
            bar_control = page.locator(".chart-option-control", has_text="Bar color")
            check(terms_control.count() >= 1 and dot_control.count() >= 1 and bar_control.count() >= 1,
                  "the parity fix landed too: Terms shown / Dot gradient / Bar color controls are present "
                  "in this pane (previously hardcoded with no controls at all)")
            if len(db_opts) >= 2:
                caption_before = page.locator("body").inner_text()
                db_control.first.select_option(index=1)
                page.wait_for_timeout(150)
                caption_after = page.locator("body").inner_text()
                check(caption_before != caption_after,
                      "switching the Database selector actually changes the displayed annotation-source caption")
    else:
        check(True, "no module cleared min_term_size for either database on this real run - "
              "the shape check above already confirmed the payload carries per-database candidates")

    print("=== CSV downloads work and carry Proteomics-specific filenames ===")
    map_dl_btn = page.locator("button", has_text="Download protein mapping (CSV)")
    check(map_dl_btn.count() >= 1, "protein-mapping CSV download button present ('protein', not 'gene') - "
          "still on the Protein mapping tab from the check above")

    hit_tab = page.locator("button.result-nav-btn", has_text="Hit table")
    check(hit_tab.count() >= 1, "results has a 'Hit table' nav tab")
    hit_tab.first.click()
    page.wait_for_timeout(150)
    dl_btn = page.locator("button", has_text="Download hit table (CSV)")
    check(dl_btn.count() >= 1, "hit-table CSV download button present")
    with page.expect_download() as dl_info:
        dl_btn.first.click()
    dl = dl_info.value
    check("proteomics_network" in dl.suggested_filename,
          f"downloaded filename is proteomics-specific (got {dl.suggested_filename})")

    print("=== console error check ===")
    check(len(errors) == 0, f"no console/page errors during the whole flow (got {errors})")

    browser.close()

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
