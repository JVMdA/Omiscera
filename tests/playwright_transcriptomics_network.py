"""Browser test for the Transcriptomics Network / Module Detection tile.

Drives the REAL flow - tile -> retriever -> resolve -> config -> run -> results
- against a live server and real network calls (mygene.info, STRING, QuickGO).
Unlike playwright_genomics_results_cards.py, which asserts only that the results
nav sits in a card, this asserts actual CONTENT: every pane renders, the hit
table carries the transcriptomics columns, recruited genes are present, module
rows are internally consistent, and the honesty caveats reach the screen.
"""
import os
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
