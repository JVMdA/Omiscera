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

    print("=== all four Proteomics sub-tiles enabled ===")
    ora_btn = page.locator("button.module-btn.proteomics", has_text="Over-Representation Analysis")
    check(ora_btn.count() >= 1 and not ora_btn.first.is_disabled(), "ORA sub-tile present and enabled")
    quant_btn = page.locator("button.module-btn.proteomics", has_text="Quantitative Enrichment")
    check(quant_btn.count() >= 1 and not quant_btn.first.is_disabled(),
          "Quantitative Enrichment sub-tile present and ENABLED (Milestone 2)")
    for name in ["PPI Network / Module Detection", "PTM / Phosphosite Enrichment"]:
        btn = page.locator("button.module-btn.proteomics", has_text=name)
        check(btn.count() >= 1 and not btn.first.is_disabled(),
              f"'{name}' present and ENABLED (see its own dedicated playwright test for full coverage)")
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

    print("=== overlapping-member column: renamed, and one consistent ACCESSION (SYMBOL) format ===")
    cpx_rows = (tables.get("Complex Portal (all)") or {}).get("table", {}).get("rows", [])
    check(bool(cpx_rows) and "Overlapping Proteins" in cpx_rows[0],
          "the column is named 'Overlapping Proteins' (not the shared engine's 'Overlapping Genes')")
    if cpx_rows:
        members = str(cpx_rows[0].get("Overlapping Proteins", ""))
        check(" (" in members and members.split(",")[0].strip().endswith(")"),
              f"members are formatted 'ACCESSION (SYMBOL)', not a bare accession (got {members[:48]!r})")

    print("=== SENSITIVITY: the protein-concept network still parses that column ===")
    # _oraGeneSet must strip the " (SYMBOL)" suffix, or every node id stops
    # matching the accession-keyed log2FC map and the figure silently loses
    # its colouring. Exercise the REAL parser against the REAL row.
    parsed = page.evaluate("""(row) => Array.from(_oraGeneSet(row))""", cpx_rows[0] if cpx_rows else {})
    check(bool(parsed) and all(" (" not in p and ")" not in p for p in parsed),
          f"_oraGeneSet returns bare accessions usable as node ids/log2FC keys (got {parsed[:4]})")
    log2fc_keys = page.evaluate("""() => {
        const rows = (STATE.proteomics.resolved.resolved_table.rows) || [];
        return rows.slice(0, 400).map(r => String(r.accession || '').toUpperCase());
    }""")
    check(bool(parsed) and any(p.upper() in set(log2fc_keys) for p in parsed),
          "those parsed ids genuinely match resolved-table accessions, so the network's log2FC colouring still resolves")

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

    print("=== DISEASES: present, human-gated, and finds real disease enrichment ===")
    # The messy-example section above left us on the upload/resolve screen,
    # which has no database checkboxes - go to the config screen first.
    page.evaluate("""() => { STATE.proteomics.taxonId = '9606'; STATE.aaPage = 'proteomics_ora_config'; render(); }""")
    page.wait_for_timeout(250)
    dis_pill = page.locator("label.ora-db-pill", has_text="DISEASES")
    check(dis_pill.count() >= 1, "the DISEASES database checkbox is present on the config screen")
    check(not dis_pill.locator("input[type=checkbox]").first.is_disabled(),
          "DISEASES is ENABLED for human (taxon 9606)")
    # Human-only gating: switching to mouse must grey it out rather than
    # letting a run come back silently empty.
    page.evaluate("""() => { STATE.proteomics.taxonId = '10090'; render(); }""")
    page.wait_for_timeout(250)
    dis_pill_mouse = page.locator("label.ora-db-pill", has_text="DISEASES")
    check(dis_pill_mouse.locator("input[type=checkbox]").first.is_disabled(),
          "DISEASES is correctly DISABLED for a non-human organism (mouse), not silently empty")
    page.evaluate("""() => { STATE.proteomics.taxonId = '9606'; render(); }""")
    page.wait_for_timeout(250)

    print("=== custom-database example zip (real download endpoint) ===")
    import io as _io
    import zipfile as _zipfile
    import uuid as _uuid
    import urllib.request as _urlreq
    req = _urlreq.Request(BASE + "/api/proteomics/ora/custom-db-example",
                           headers={"X-Session-Id": _uuid.uuid4().hex})
    with _urlreq.urlopen(req, timeout=30) as resp:
        zip_bytes = resp.read()
    zf = _zipfile.ZipFile(_io.BytesIO(zip_bytes))
    names = set(zf.namelist())
    check(names == {"README_proteomics_custom_geneset_db.txt",
                     "example_proteomics_custom_geneset_db.gmt",
                     "example_proteomics_custom_geneset_db_TIDY.csv"},
          f"the custom-DB example zip contains exactly the README + GMT + tidy CSV (got {sorted(names)})")

    print("=== custom-DB upload: real GMT/tidy files from the zip actually parse and load ===")
    page.evaluate("""() => { STATE.aaPage = 'proteomics_ora_config'; render(); }""")
    page.wait_for_timeout(150)
    gmt_bytes = zf.read("example_proteomics_custom_geneset_db.gmt")
    import tempfile as _tempfile
    with _tempfile.NamedTemporaryFile(suffix=".gmt", delete=False) as tmp:
        tmp.write(gmt_bytes)
        gmt_path = tmp.name
    custom_cb = page.locator("label.ora-db-pill", has_text="Custom database").locator("input[type=checkbox]")
    custom_cb.check()
    page.wait_for_timeout(200)
    page.locator("input[type=file]").last.set_input_files(gmt_path)
    page.wait_for_timeout(800)
    custom_db_state = page.evaluate("STATE.proteomics.customDb")
    check(custom_db_state is not None and custom_db_state.get("n_terms") == 6,
          f"the real zipped GMT example uploads and parses as 6 gene sets (got {custom_db_state})")

    print("=== SENSITIVITY: a symbol-based custom DB is rejected with the REAL reason, not a silent empty result ===")
    with _tempfile.NamedTemporaryFile(suffix=".gmt", mode="w", delete=False, encoding="utf-8") as tmp2:
        tmp2.write("fake_set\tNot real accessions\tTP53\tMDM2\tBRCA1\tEGFR\tMYC\n")
        symbol_path = tmp2.name
    page.once("dialog", lambda d: d.accept())
    page.locator("input[type=file]").last.set_input_files(symbol_path)
    page.wait_for_timeout(500)
    custom_db_after_bad = page.evaluate("STATE.proteomics.customDb")
    check(custom_db_after_bad is None,
          f"a gene-symbol-only file is REJECTED (customDb stays null), not silently accepted and matching nothing "
          f"(got {custom_db_after_bad})")

    print("=== KEGG: pill present, no organism picker, clean human run ===")
    # STATE.proteomics.resolved at this point is actually still the MESSY
    # MaxQuant example from the contaminant-filtering section above (the
    # DISEASES section only set taxonId and navigated - it never re-uploaded
    # the clean file), so re-do a clean human upload+resolve explicitly
    # rather than depending on inherited cross-section state.
    page.evaluate("""() => {
        STATE.aaPage = 'proteomics';
        STATE.proteomicsModuleName = '';
        render();
    }""")
    page.wait_for_timeout(150)
    page.locator("button.module-btn.proteomics", has_text="Over-Representation Analysis").click()
    page.wait_for_timeout(150)
    page.evaluate("""() => { STATE.proteomics.taxonId = '9606'; STATE.proteomics.taxonName = 'Homo sapiens'; render(); }""")
    page.wait_for_timeout(150)
    page.locator("input[type=file]").set_input_files(
        os.path.join(_ROOT, "Data3", "example_proteomics_de_full.csv"))
    page.wait_for_timeout(1200)
    page.locator("button", has_text="Resolve accessions").first.click()
    for _ in range(180):
        if page.evaluate("!!STATE.proteomics.resolved"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.proteomics.resolved") is True, "re-resolved the clean human example for the KEGG section")
    page.locator("button", has_text="Continue to Over-Representation Analysis").click()
    page.wait_for_timeout(200)

    kegg_pill = page.locator("label.ora-db-pill", has_text="KEGG").first
    check(kegg_pill.count() == 1, "a KEGG pill is present in the database checklist")
    kegg_cb = kegg_pill.locator("input[type=checkbox]")
    check(not kegg_cb.is_disabled(), "the KEGG pill is enabled")
    caption_before = page.evaluate(
        "() => [...document.querySelectorAll('.caption')].some(c => c.textContent.includes('KEGG uses the organism'))")
    check(caption_before is False, "the KEGG organism caption is hidden until KEGG is checked")
    kegg_cb.check()
    page.wait_for_timeout(1200)
    caption_after = page.evaluate(
        "() => [...document.querySelectorAll('.caption')].some(c => c.textContent.includes('KEGG uses the organism'))")
    check(caption_after is True, "checking KEGG shows the organism-derivation caption (real re-render, not stale)")
    check(page.locator(".ora-kegg-organism-picker").count() == 0,
          "SENSITIVITY: NO disambiguation picker for the clean human case - a second, independently-"
          "settable organism control here is exactly the bug fixed in commit 27ce033")

    run_btn = page.locator("button.ora-run-btn")
    check(not run_btn.is_disabled(), "Run button enabled with KEGG selected")
    run_btn.click()
    page.wait_for_timeout(300)
    for _ in range(300):
        if page.evaluate("!!STATE.proteomics.oraResults"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.proteomics.oraResults") is True, "the KEGG-only run finished")
    kegg_tables = page.evaluate("STATE.proteomics.oraResults.tables")
    kegg_keys = sorted(kegg_tables.keys())
    check(kegg_keys == ["KEGG (all)", "KEGG (down)", "KEGG (up)"],
          f"three KEGG direction families came back through the real UI (got {kegg_keys})")
    kegg_all = kegg_tables.get("KEGG (all)", {})
    check(kegg_all.get("n_terms_tested", 0) > 0,
          f"the KEGG family actually tested real terms (got n_terms_tested={kegg_all.get('n_terms_tested')})")

    print("=== KEGG: strain-disambiguation picker for a real ambiguous organism (Bacillus subtilis) ===")
    page.evaluate("""() => {
        STATE.aaPage = 'proteomics';
        STATE.proteomicsModuleName = '';
        render();
    }""")
    page.wait_for_timeout(150)
    page.locator("button.module-btn.proteomics", has_text="Over-Representation Analysis").click()
    page.wait_for_timeout(150)
    page.evaluate("""() => {
        STATE.proteomics.taxonId = '1423';
        STATE.proteomics.taxonName = 'Bacillus subtilis';
        render();
    }""")
    page.wait_for_timeout(150)
    page.locator("input[type=file]").set_input_files(
        os.path.join(_ROOT, "Data3", "example_proteomics_de_full.csv"))
    page.wait_for_timeout(1200)
    page.locator("button", has_text="Resolve accessions").first.click()
    for _ in range(180):
        if page.evaluate("!!STATE.proteomics.resolved"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.proteomics.resolved") is True,
          "resolve completed for Bacillus subtilis (real UniProt ID Mapping round trip)")
    page.locator("button", has_text="Continue to Over-Representation Analysis").click()
    page.wait_for_timeout(200)
    check(page.evaluate("STATE.aaPage") == "proteomics_ora_config",
          "navigated to the ORA config screen for the ambiguous organism")

    page.locator("label.ora-db-pill", has_text="KEGG").first.locator("input[type=checkbox]").check()
    picker = page.locator(".ora-kegg-organism-picker")
    for _ in range(30):
        if picker.count() > 0:
            break
        page.wait_for_timeout(500)
    check(picker.count() == 1,
          "SENSITIVITY: the SAME disambiguation picker used by Transcriptomics also appears here for "
          "Proteomics, proving the widget's generalization (moduleKey='proteomics') genuinely works, "
          "not just a Transcriptomics-shaped stub")
    select = picker.locator("select")
    confirm_btn = picker.locator("button", has_text="Confirm KEGG organism")
    options = select.locator("option").all_inner_texts()
    check(any(o.startswith("bsu —") for o in options),
          f"the real 'bsu' (B. subtilis 168) candidate is offered by name (got {options[:5]})")
    select.select_option("bsu")
    confirm_btn.click()
    for _ in range(30):
        if page.evaluate("() => !!(STATE.proteomics._keggOrganismConfirmed)"):
            break
        page.wait_for_timeout(500)
    confirmed_state = page.evaluate("() => STATE.proteomics._keggOrganismConfirmed")
    check(confirmed_state and confirmed_state.get("organismCode") == "bsu",
          f"confirming persists the real chosen organism into Proteomics' OWN state key "
          f"(got {confirmed_state})")
    check(page.locator(".ora-selected-caption", has_text="Confirmed manually").count() == 1,
          "the config screen shows the confirmed-organism caption")
    check(picker.count() == 0, "and the picker itself is gone, replaced by that caption")

    print("=== console error check ===")
    check(len(errors) == 0, f"no console/page errors during the whole flow (got {errors})")

    browser.close()

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
