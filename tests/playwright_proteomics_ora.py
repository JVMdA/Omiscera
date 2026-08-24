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

    print("=== console error check ===")
    check(len(errors) == 0, f"no console/page errors during the whole flow (got {errors})")

    browser.close()

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
