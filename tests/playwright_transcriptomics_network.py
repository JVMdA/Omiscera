"""Browser test for the Transcriptomics Network / Module Detection tile.

Drives the REAL flow - tile -> retriever -> resolve -> config -> run -> results
- against a live server and real network calls (mygene.info, STRING, QuickGO).
Unlike playwright_genomics_results_cards.py, which asserts only that the results
nav sits in a card, this asserts actual CONTENT: every pane renders, the hit
table carries the transcriptomics columns, recruited genes are present, module
rows are internally consistent, and the honesty caveats reach the screen.
"""
import os
import re
import sys

from playwright.sync_api import sync_playwright

BASE = "http://127.0.0.1:8009"
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) \
    if os.path.basename(os.path.dirname(os.path.abspath(__file__))) == "tests" \
    else r"c:\Users\cliente\OneDrive\Documentos\Doutorado_brasil\Resultados\app-metabolismo\Código\Omiscera_v06"
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
    page.on("pageerror", lambda e: errors.append(str(e)))
    page.on("console", lambda m: errors.append("console: " + m.text) if m.type == "error" else None)
    page.goto(BASE, wait_until="networkidle")
    page.wait_for_timeout(300)

    check(page.evaluate("typeof renderTranscriptomicsNetworkConfigScreen") == "function",
          "the network config screen is defined (app.js parsed with no top-level error)")

    print("=== the fourth tile is present and enabled ===")
    page.evaluate("""() => { STATE.view='annotation_analytics'; STATE.aaPage='transcriptomics'; render(); }""")
    page.wait_for_timeout(200)
    tile = page.locator("button.module-btn.transcriptomics", has_text="Network / Module Detection")
    check(tile.count() == 1, "the 'Network / Module Detection' tile is present")
    check(not tile.is_disabled(), "the tile is ENABLED (no longer the disabled 'Enrichment (coming soon)')")
    labels = page.evaluate("[...document.querySelectorAll('button.module-btn.transcriptomics')].map(b=>b.textContent)")
    check(not any("coming soon" in t for t in labels),
          f"no Transcriptomics tile is still marked 'coming soon' (got {labels})")
    tile.click()
    page.wait_for_timeout(200)

    print("=== retriever: the organism selector searches STRING for this tile ===")
    # The organism selector only renders once an identifier type is chosen.
    page.locator("button", has_text="Gene Symbol").first.click()
    page.wait_for_timeout(250)
    ph = page.locator("input.text-input[placeholder]").first.get_attribute("placeholder")
    check("STRING" in (ph or ""),
          f"the organism box says it searches STRING, not the whole UniProt taxonomy (got {ph!r})")
    page.evaluate("""() => { STATE.transcriptomics.taxonId='9606'; STATE.transcriptomics.taxonName='Homo sapiens'; render(); }""")
    page.wait_for_timeout(150)
    page.locator("input[type=file]").set_input_files(
        os.path.join(_ROOT, "Data3", "example_transcriptomics_de_full.csv"))
    page.wait_for_timeout(1500)
    page.locator("button", has_text="Resolve identifiers").click()
    for _ in range(180):
        if page.evaluate("!!STATE.transcriptomics.resolved"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.transcriptomics.resolved") is True, "resolve completed (real mygene.info)")

    print("=== continue -> config -> run ===")
    cont = page.locator("button", has_text="Continue to Network / Module Detection")
    check(cont.count() == 1, "the Continue button routes to this module")
    cont.click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "transcriptomics_network_config", "navigated to the config screen")

    # Keep the run fast: fewest permutations offered.
    page.evaluate("""() => { STATE.transcriptomics._networkDraft.nPermutations = 200; render(); }""")
    page.wait_for_timeout(200)
    run = page.locator("button.ora-run-btn")
    check(run.count() == 1 and not run.is_disabled(), "the Run button is enabled")
    run.click()
    for _ in range(600):
        if page.evaluate("!!STATE.transcriptomics.networkResults"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.transcriptomics.networkResults") is True, "the network job completed")
    page.wait_for_timeout(500)

    print("=== results: panes, card, content ===")
    rails = page.evaluate("[...document.querySelectorAll('.result-nav-btn')].map(b=>b.textContent)")
    for want in ("Overview", "Gene mapping", "Hit table", "Figures", "Parameters used"):
        check(want in rails, f"'{want}' pane is present (rails: {rails})")
    check(page.evaluate("""() => {
        const nav = document.querySelector('.result-nav-layout');
        return !!(nav && nav.closest('.ora-db-card'));
    }"""), "the results nav sits inside the standard '.ora-db-card' box")

    stats = page.evaluate("""() => {
        const p = STATE.transcriptomics.networkResults.propagation;
        return { hits: p.n_hits, recruited: p.n_recruited, seedHits: p.n_hits_seed,
                 cols: p.table.columns, edges: (p.edges||[]).length };
    }""")
    check(stats["hits"] > 0, f"the run produced significant proteins ({stats['hits']})")
    check(stats["recruited"] > 0,
          f"including RECRUITED proteins never uploaded ({stats['recruited']}) - what this module adds")
    check("log2FC" in stats["cols"] and "Mapped via" in stats["cols"] and "UniProt ID" not in stats["cols"],
          f"the hit table is transcriptomics-shaped, not the Genomics UniProt one ({stats['cols']})")

    caveats = page.evaluate("[...document.querySelectorAll('.alert-warning')].map(e=>e.textContent)")
    check(any("seeds" in c for c in caveats),
          "the seed-circularity caveat is shown on screen, not buried in the payload")

    print("=== evaluation Part 3c: the on-screen caveat's counts stay honest ===")
    seed_caveat = next((c for c in caveats if "supplied as seeds" in c), None)
    check(seed_caveat is not None, "the seed-circularity caveat text is findable for parsing")
    if seed_caveat:
        m = re.search(r"(\d+) of the (\d+) significant proteins are genes you supplied as seeds", seed_caveat)
        m2 = re.search(r"The (\d+) RECRUITED protein", seed_caveat)
        check(bool(m and m2), f"the caveat text still matches the expected pattern (got: {seed_caveat!r})")
        if m and m2:
            check(int(m.group(1)) == stats["seedHits"] and int(m.group(2)) == stats["hits"]
                  and int(m2.group(1)) == stats["recruited"],
                  f"the caveat's own numbers (seed={m.group(1)}, total={m.group(2)}, recruited={m2.group(1)}) "
                  f"match the payload's n_hits_seed/n_hits/n_recruited exactly (payload: "
                  f"{stats['seedHits']}/{stats['hits']}/{stats['recruited']}) - the Part 3b toggle only "
                  f"filters the Hit table's DISPLAY, it never touches these counts")

    print("=== evaluation Part 3a: the shared FDR-underflow formatter is wired in, on real data ===")
    fmt = page.evaluate("""() => ({
        zero: _fmtSci(0, 3), empty: _fmtSci("", 3), nullVal: _fmtSci(null, 3),
        normal: _fmtSci(0.001234, 3), fixedEmpty: _fmtFixed("", 3), fixedNormal: _fmtFixed(2.71828, 3),
    })""")
    check(fmt["zero"] == "< 1e-300", f"_fmtSci(0) renders '< 1e-300', not the pre-fix '0.000e+0' (got {fmt['zero']!r})")
    check(fmt["empty"] == "" and fmt["nullVal"] == "",
          f"a MISSING value (a recruited gene's blank padj) stays blank, not '< 1e-300' (got {fmt['empty']!r}, {fmt['nullVal']!r})")
    check(fmt["normal"] == "1.234e-3", f"a normal p-value still formats in scientific notation (got {fmt['normal']!r})")
    check(fmt["fixedEmpty"] == "" and fmt["fixedNormal"] == "2.718",
          f"_fmtFixed also blanks missing values and rounds normal ones (got {fmt['fixedEmpty']!r}, {fmt['fixedNormal']!r})")

    real_underflow = page.evaluate("""() => {
        const rows = STATE.transcriptomics.networkResults.propagation.table.rows || [];
        return rows.some(r => r["FDR"] === 0);
    }""")
    check(real_underflow, "this real run genuinely produced an FDR of exactly 0.0 (not simulated) - "
                          "the case the fix above exists for")
    page.locator(".result-nav-btn", has_text="Hit table").first.click()
    page.wait_for_timeout(300)
    cell_texts = page.evaluate("[...document.querySelectorAll('.data-table td')].map(td => td.textContent)")
    check("< 1e-300" in cell_texts, "and the on-screen Hit table actually renders it as '< 1e-300'")
    check(not any(t.startswith("0.000e") for t in cell_texts),
          "and no cell anywhere in the Hit table shows the broken '0.000e+0'")

    csv = page.evaluate("""async () => {
        const r = await fetch('/api/transcriptomics/network/download?which=hits',
                              { headers: { 'X-Session-Id': getSessionId() } });
        return await r.text();
    }""")
    check("< 1e-300" not in csv,
          "the CSV download keeps full numeric precision - it never contains the display-only '< 1e-300' text")

    print("=== evaluation Part 3b: the 'recruited only' Hit table toggle ===")
    toggle_state = page.evaluate("""() => {
        const cb = [...document.querySelectorAll('input[type=checkbox]')]
            .find(c => c.parentElement.textContent.includes('RECRUITED genes only'));
        return cb ? cb.checked : null;
    }""")
    check(toggle_state is True,
          "with real recruited genes present, the toggle defaults to ON - the discovery is shown, not buried")
    def _hit_pane_row_count():
        return page.evaluate("""() => {
            const sections = [...document.querySelectorAll('.result-pane')];
            const hitSection = sections.find(
                s => s.querySelector('.result-pane-title') && s.querySelector('.result-pane-title').textContent === 'Hit table');
            return hitSection ? hitSection.querySelectorAll('.data-table tbody tr').length : -1;
        }""")

    n_recruited = page.evaluate("STATE.transcriptomics.networkResults.propagation.n_recruited")
    rows_checked = _hit_pane_row_count()
    check(rows_checked == n_recruited,
          f"checked, the table shows EXACTLY the {n_recruited} recruited rows, not the full {stats['hits']}-row "
          f"significant set (got {rows_checked} rows)")
    checked_seed_cells = page.evaluate("""() => {
        const sections = [...document.querySelectorAll('.result-pane')];
        const hitSection = sections.find(
            s => s.querySelector('.result-pane-title') && s.querySelector('.result-pane-title').textContent === 'Hit table');
        const rows = [...hitSection.querySelectorAll('.data-table tbody tr')];
        return rows.some(tr => tr.textContent.includes('true'));
    }""")
    check(not checked_seed_cells, "no seed row leaks into the recruited-only view")

    page.evaluate("""() => {
        const cb = [...document.querySelectorAll('input[type=checkbox]')]
            .find(c => c.parentElement.textContent.includes('RECRUITED genes only'));
        cb.checked = false;
        cb.dispatchEvent(new Event('change'));
    }""")
    page.wait_for_timeout(200)
    rows_unchecked = _hit_pane_row_count()
    check(rows_unchecked > rows_checked,
          f"unchecked, the SAME table now shows the full network ({rows_unchecked} rows, up from {rows_checked}) - "
          f"the toggle only filters the view, it does not refetch or discard data")

    dl_csv_after_toggle = page.evaluate("""async () => {
        const r = await fetch('/api/transcriptomics/network/download?which=hits',
                              { headers: { 'X-Session-Id': getSessionId() } });
        const t = await r.text();
        return t.split('\\n').length - 1;
    }""")
    check(dl_csv_after_toggle > rows_unchecked,
          f"the CSV download still carries EVERY protein in the network ({dl_csv_after_toggle} rows), regardless "
          f"of the on-screen toggle's state")

    print("=== figures render ===")
    page.locator(".result-nav-btn", has_text="Figures").first.click()
    page.wait_for_timeout(300)
    page.locator("button.chart-toggle-btn", has_text="Show network figure").first.click()
    page.wait_for_timeout(2500)
    net = page.evaluate("""() => {
        const d = document.querySelector('#chart-propnet-transcriptomics_network-score');
        return d && d.data ? d.data.length : 0;
    }""")
    check(net >= 2, f"the network figure rendered edge+node traces ({net})")
    page.locator("button.chart-toggle-btn", has_text="Show null distribution").first.click()
    page.wait_for_timeout(1500)
    hist = page.evaluate("""() => {
        const d = document.querySelector('#chart-histogram-transcriptomics_network-null');
        return d && d.data ? d.data.length : 0;
    }""")
    check(hist >= 1, f"the null-distribution figure rendered ({hist})")

    mods = page.evaluate("""() => {
        const m = STATE.transcriptomics.networkResults.modules;
        if (!m) return null;
        const rows = m.table.rows || [];
        return { n: rows.length,
                 consistent: rows.every(r => (r['Seeds'] + r['Recruited']) === r['Size']) };
    }""")
    if mods:
        check(mods["n"] > 0, f"modules were detected ({mods['n']})")
        check(mods["consistent"] is True,
              "every module row is internally consistent: Seeds + Recruited == Size")
        check("Modules" in rails, "the Modules pane is present")

    check(len(errors) == 0, f"no console/page errors during the whole flow (got {errors[:3]})")
    browser.close()

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
