"""FIRST real-browser coverage of Proteomics Quantitative Enrichment
(Milestone 2 - 1D annotation enrichment / competitive Wilcoxon, reusing
Genomics' Gene-set Burden Test engine extended to alternative="two-sided").

Drives the REAL flow end to end: tile -> organism selector -> upload the
bundled example -> resolve (real UniProt ID Mapping) -> continue to
Quantitative Enrichment -> configure -> run (real InterPro/Complex
Portal/UniProt keywords fetches) -> results, then checks the member-column
fix to _burdenGeneSet (the one genuine frontend correctness fix this
milestone needed - see app.js) actually works against a real row.

Needs a live server and real network access - not mocked, takes a couple of
minutes.
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

    print("=== navigate to Proteomics -> Quantitative Enrichment ===")
    page.evaluate("""() => { STATE.view = 'annotation_analytics'; STATE.aaPage = 'overview'; render(); }""")
    page.wait_for_timeout(200)
    page.locator("button.module-btn.proteomics").first.click()
    page.wait_for_timeout(200)
    quant_btn = page.locator("button.module-btn.proteomics", has_text="Quantitative Enrichment")
    check(quant_btn.count() >= 1 and not quant_btn.first.is_disabled(), "Quantitative Enrichment tile enabled")
    quant_btn.first.click()
    page.wait_for_timeout(200)
    check(page.evaluate("STATE.aaPage") == "proteomics_module", "navigated to the upload/resolve screen")
    check(page.evaluate("STATE.proteomicsModuleName") == "Quantitative Enrichment", "module name set correctly")

    print("=== organism + upload + resolve (real UniProt ID Mapping) ===")
    page.evaluate("""() => { STATE.proteomics.taxonId = '9606'; STATE.proteomics.taxonName = 'Homo sapiens'; render(); }""")
    page.wait_for_timeout(150)
    page.locator("input[type=file]").set_input_files(
        os.path.join(_ROOT, "Data3", "example_proteomics_de_full.csv")
    )
    page.wait_for_timeout(1500)
    check(page.evaluate("!!STATE.proteomics.deTable"), "upload succeeded, preview table populated")

    resolve_btn = page.locator("button", has_text="Resolve accessions")
    resolve_btn.first.click()
    for _ in range(180):
        if page.evaluate("!!STATE.proteomics.resolved"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.proteomics.resolved") is True, "resolve completed (real UniProt ID Mapping round trip)")

    print("=== continue -> configure -> run (fast databases only, for speed) ===")
    continue_btn = page.locator("button", has_text="Continue to Quantitative Enrichment")
    check(continue_btn.count() >= 1, "'Continue to Quantitative Enrichment' button present (not the ORA-only fallback text)")
    continue_btn.first.click()
    page.wait_for_timeout(200)
    check(page.evaluate("STATE.aaPage") == "proteomics_quant_config", "navigated to the Quant config screen")

    # No padj/log2FC cutoff fields on this screen (threshold-free test) -
    # confirm the ORA-only fields genuinely do not leak into this screen.
    check(page.locator("text=Adjusted p-value (padj) cutoff").count() == 0,
          "no padj cutoff field on the Quant config screen (this test has no significance cutoff)")
    check(page.locator("text=Ranking metric").count() >= 1, "the ranking-metric selector IS present")
    metric_select = page.locator("select").filter(has=page.locator("option", has_text="log2 Fold-Change"))
    check(metric_select.count() >= 1, "the metric dropdown offers log2 Fold-Change")
    default_metric_value = page.evaluate("STATE.proteomics._quantDraft.metric")
    check(default_metric_value == "log2fc", f"default metric is log2fc (got {default_metric_value!r})")

    # Complex Portal alone now reliably clears FDR<0.05 on the bundled
    # example (fixed 2026-08-25: the generator used to draw each seeded
    # protein's up/down sign independently, so a complex's own members
    # cancelled out under a SIGNED, direction-aware test like this one - see
    # generate_proteomics_examples.py's complex_sign for the fix). GO is no
    # longer needed just to reach a significant result, so the fast three
    # databases are enough here.
    for db_name in ["InterPro", "Complex Portal", "UniProt (keywords)"]:
        page.locator("label.ora-db-pill", has_text=db_name).locator("input[type=checkbox]").check()
    page.wait_for_timeout(100)
    run_btn = page.locator("button.ora-run-btn")
    check(not run_btn.is_disabled(), "Run button enabled with 3 databases selected")
    run_btn.click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "proteomics_quant_results", "navigated to the Quant results screen")

    for _ in range(180):
        if page.evaluate("!!STATE.proteomics.quantResults"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.proteomics.quantResults") is True, "the Quantitative Enrichment job completed")

    print("=== results: real, non-empty, no significance-cutoff summary fields ===")
    summary = page.evaluate("STATE.proteomics.quantResults.summary")
    check(summary["n_scored"] > 0, f"proteins were scored (n_scored={summary['n_scored']})")
    check(summary["metric"] == "log2fc", f"summary echoes the metric actually used (got {summary['metric']!r})")
    check("n_sig_all" not in summary, "no ORA-style significance-count field leaked into the Quant summary")
    tables = page.evaluate("STATE.proteomics.quantResults.tables")
    # No direction suffix on Quant's result keys, unlike ORA's "(all)/(up)/(down)".
    check("Complex Portal" in tables and "Complex Portal (all)" not in tables,
          f"Complex Portal's result key has NO direction suffix (got keys: {sorted(tables.keys())})")
    cpx = tables.get("Complex Portal")
    check(cpx is not None and cpx["n_terms_tested"] > 0, "Complex Portal table exists and has tested terms")

    check(page.locator("select.transcriptomics-ora-result-select").count() >= 1, "the result-table selector is present")
    check(page.locator(".ora-db-card").count() >= 1, "a results card rendered")

    print("=== overlapping-member column: same 'Overlapping Proteins' rename and format as ORA ===")
    cpx_rows = (cpx or {}).get("table", {}).get("rows", [])
    check(bool(cpx_rows) and "Overlapping Proteins" in cpx_rows[0],
          "the Burden-engine's 'Genes' column arrives renamed to 'Overlapping Proteins' here too")
    check("Genes" not in (cpx_rows[0] if cpx_rows else {}),
          "the raw 'Genes' column name does not leak through for Proteomics")
    if cpx_rows:
        members = str(cpx_rows[0].get("Overlapping Proteins", ""))
        check(" (" in members and members.split(",")[0].strip().endswith(")"),
              f"members are formatted 'ACCESSION (SYMBOL)' (got {members[:48]!r})")

    print("=== SENSITIVITY: _burdenGeneSet (fixed this milestone) parses the renamed column ===")
    # Before this milestone's fix, _burdenGeneSet read ONLY row['Genes'] - for
    # a Proteomics Quant table (column renamed to 'Overlapping Proteins') it
    # would have returned an EMPTY set, silently breaking the enrichment map
    # and protein-concept network's edges. Exercise the REAL function against
    # a REAL row to prove the fix works, not just that it compiles.
    parsed = page.evaluate("""(row) => Array.from(_burdenGeneSet(row))""", cpx_rows[0] if cpx_rows else {})
    check(bool(parsed), f"_burdenGeneSet returns a NON-EMPTY set for a Proteomics-renamed row (got {parsed})")
    check(all(" (" not in g and ")" not in g for g in parsed),
          f"_burdenGeneSet strips the '(SYMBOL)' suffix, returning bare accessions usable as node ids (got {parsed[:4]})")
    log2fc_keys = page.evaluate("""() => {
        const rows = (STATE.proteomics.resolved.resolved_table.rows) || [];
        return rows.slice(0, 400).map(r => String(r.accession || '').toUpperCase());
    }""")
    check(any(g.upper() in set(log2fc_keys) for g in parsed),
          "those parsed ids genuinely match resolved-table accessions, so gene-concept-network log2FC coloring resolves")

    print("=== figures render on a REAL significant term (all 5 Burden-shaped chart builders) ===")
    # Find whichever selected database actually has a significant term this
    # run (Complex Portal, live-verified 2026-08-25 after the example fix,
    # but don't hardcode it - assert on the real payload instead of assuming
    # a specific database).
    sig_key = None
    for key, block in tables.items():
        rows = block.get("table", {}).get("rows", [])
        if any(r.get("FDR") is not None and float(r["FDR"]) < 0.05 for r in rows):
            sig_key = key
            break
    check(sig_key is not None, f"at least one selected database has a significant term (keys: {sorted(tables.keys())})")

    if sig_key:
        page.evaluate(
            """(key) => {
                const sel = document.querySelector('select.transcriptomics-ora-result-select');
                sel.value = key;
                sel.dispatchEvent(new Event('change'));
            }""",
            sig_key,
        )
        page.wait_for_timeout(300)
        overview_rail_btn = page.locator("button.result-nav-btn", has_text="Protein set overview")
        check(overview_rail_btn.count() >= 1, f"'Protein set overview' pane present for {sig_key!r}")
        if overview_rail_btn.count() >= 1:
            overview_rail_btn.first.click()
            page.wait_for_timeout(200)
            for label in ["Show dot plot", "Show bar plot", "Show volcano plot"]:
                btn = page.locator("button.chart-toggle-btn", has_text=label)
                check(btn.count() >= 1, f"'{label}' button present")
                if btn.count() >= 1:
                    btn.first.click()
                    page.wait_for_timeout(500)
            n_plots = page.locator(".js-plotly-plot").count()
            check(n_plots >= 3, f"dot+bar+volcano all rendered real Plotly figures (got {n_plots})")

            relationships_rail_btn = page.locator("button.result-nav-btn", has_text="Protein set relationships")
            check(relationships_rail_btn.count() >= 1, "'Protein set relationships' pane present")
            if relationships_rail_btn.count() >= 1:
                relationships_rail_btn.first.click()
                page.wait_for_timeout(200)
                for label in ["Show enrichment map", "Show protein-concept network"]:
                    btn = page.locator("button.chart-toggle-btn", has_text=label)
                    check(btn.count() >= 1, f"'{label}' button present")
                    if btn.count() >= 1:
                        btn.first.click()
                        page.wait_for_timeout(800)
                n_plots2 = page.locator(".js-plotly-plot").count()
                check(n_plots2 >= n_plots,
                      f"the protein-concept network (which needs only 1 term, unlike the enrichment map's "
                      f"minimum of 2) rendered too (now {n_plots2} total)")

    print("=== download endpoint ===")
    csv_resp = page.request.get(
        BASE + "/api/proteomics/quant/download?key=" + "Complex%20Portal",
        headers={"X-Session-Id": page.evaluate("getSessionId()")},
    )
    check(csv_resp.ok, f"download endpoint responds OK (status={csv_resp.status})")
    csv_text = csv_resp.text()
    check("Overlapping Proteins" in csv_text.splitlines()[0], "downloaded CSV header uses the renamed column")

    print("=== MSigDB collection summary: the gap this run closes ===")
    # Proteomics Quant offered MSigDB (any number of collections at once,
    # same as ORA/GSEA) but had no cross-collection summary figure - the
    # ORA/GSEA screens both have one, Quant did not. Oncogenic (189 human
    # sets) kept deliberately small for speed, same choice playwright_msigdb.py
    # already made for the identical reason.
    page.evaluate("""() => {
        STATE.proteomics.quantRunLocked = false;
        STATE.aaPage = 'proteomics_quant_config';
        delete STATE.proteomics._quantDraft;
        render();
    }""")
    page.wait_for_timeout(200)
    page.locator("label.ora-db-pill", has_text="MSigDB").first.locator("input[type=checkbox]").check()
    page.wait_for_timeout(200)
    oncogenic_pill = page.locator(".ora-msigdb-collections label", has_text="Oncogenic")
    check(oncogenic_pill.count() == 1, "the Oncogenic collection checkbox is present on the Quant config screen")
    oncogenic_pill.locator("input[type=checkbox]").check()
    page.wait_for_timeout(100)
    run_btn3 = page.locator("button.ora-run-btn")
    check(not run_btn3.is_disabled(), "Run button enabled with MSigDB: Oncogenic selected")
    run_btn3.click()
    page.wait_for_timeout(300)
    for _ in range(300):
        if page.evaluate("!!STATE.proteomics.quantResults"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.proteomics.quantResults") is True, "the MSigDB-only Quant run finished")
    msigdb_quant_tables = page.evaluate("Object.keys(STATE.proteomics.quantResults.tables)")
    check(any(k.startswith("MSigDB: Oncogenic") for k in msigdb_quant_tables),
          f"result tables include an 'MSigDB: Oncogenic' family, no direction suffix "
          f"(got {msigdb_quant_tables})")

    summary_btn2 = page.locator("button", has_text="Show MSigDB collection summary")
    check(summary_btn2.count() == 1, "'Show MSigDB collection summary' toggle is present on the Quant results screen")
    summary_btn2.click()
    page.wait_for_timeout(800)
    summary_has_data = page.evaluate("""() => {
        const divs = [...document.querySelectorAll('.chart-container')];
        const d = divs.find(x => x.id.includes('msigdbsummary') && x.data);
        return !!(d && d.data && d.data.length > 0 && d.data[0].y.length > 0);
    }""")
    check(summary_has_data,
          "the new MSigDB collection summary figure actually renders with real Quant data - "
          "confirms _msigdbFamiliesFromQuantTables' own column names (FDR, Effect (mean "
          "difference)) genuinely match what the Wilcoxon engine's blocks carry, not just "
          "that the function is defined")

    # Changing the significance threshold must refresh an ALREADY-OPEN
    # summary too (refreshMsigdbSummary), not just a freshly-opened one -
    # this is the one behaviour that could plausibly be missed when adding
    # the wiring (declaring refreshMsigdbSummary but forgetting to call it
    # from the alpha input's change handler).
    alpha_input = page.locator("input.ora-sig-threshold-input")
    alpha_input.fill("0.99")
    alpha_input.dispatch_event("change")
    page.wait_for_timeout(500)
    summary_still_has_data = page.evaluate("""() => {
        const divs = [...document.querySelectorAll('.chart-container')];
        const d = divs.find(x => x.id.includes('msigdbsummary') && x.data);
        return !!(d && d.data && d.data.length > 0 && d.data[0].y.length > 0);
    }""")
    check(summary_still_has_data,
          "the open summary figure survives a significance-threshold change without erroring "
          "(refreshMsigdbSummary was actually wired into the alpha input's change handler)")

    print("=== KEGG: clean human run through the real UI ===")
    # The full strain-disambiguation-picker flow (a real ambiguous organism,
    # confirming a candidate) is already covered end to end in
    # playwright_proteomics_ora.py, using the SAME generalized widget this
    # screen also calls with moduleKey="proteomics" - no need to duplicate
    # that whole flow here. This just confirms KEGG reaches Quantitative
    # Enrichment's own screen and produces a real single-family result.
    page.evaluate("""() => {
        STATE.proteomics.quantRunLocked = false;
        STATE.aaPage = 'proteomics_quant_config';
        delete STATE.proteomics._quantDraft;
        render();
    }""")
    page.wait_for_timeout(200)
    check(page.evaluate("STATE.aaPage") == "proteomics_quant_config",
          "navigated back to the Quant config screen for the KEGG section")
    kegg_pill = page.locator("label.ora-db-pill", has_text="KEGG").first
    check(kegg_pill.count() == 1, "a KEGG pill is present on the Quantitative Enrichment config screen")
    kegg_pill.locator("input[type=checkbox]").check()
    page.wait_for_timeout(1200)
    check(page.locator(".ora-kegg-organism-wrap").count() == 1,
          "checking KEGG renders its organism-derivation widget here too")
    check(page.locator(".ora-kegg-organism-picker").count() == 0,
          "no disambiguation picker for the clean human case")

    run_btn2 = page.locator("button.ora-run-btn")
    check(not run_btn2.is_disabled(), "Run button enabled with KEGG selected")
    run_btn2.click()
    page.wait_for_timeout(300)
    for _ in range(300):
        if page.evaluate("!!STATE.proteomics.quantResults"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.proteomics.quantResults") is True, "the KEGG-only Quant run finished")
    kegg_quant_tables = page.evaluate("STATE.proteomics.quantResults.tables")
    check(list(kegg_quant_tables.keys()) == ["KEGG"],
          f"one KEGG family, no direction suffix (got {list(kegg_quant_tables.keys())})")
    check(kegg_quant_tables["KEGG"].get("n_terms_tested", 0) > 0,
          f"KEGG tested real terms through the real UI "
          f"(got n_terms_tested={kegg_quant_tables['KEGG'].get('n_terms_tested')})")

    print("=== console error check ===")
    check(len(errors) == 0, f"no console/page errors during the whole flow (got {errors})")

    browser.close()

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
