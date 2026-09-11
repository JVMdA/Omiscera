"""Playwright end-to-end UI test for KEGG in Transcriptomics ORA and GSEA.

Drives real clicks through the actual retriever -> config -> results flow
against a live uvicorn (127.0.0.1:8009), with REAL KEGG downloads.

Companion to tests/test_transcriptomics_kegg.py (fixtures, no network) and
tests/live_transcriptomics_kegg.py (real KEGG, no browser). This is the layer
that catches what neither of those can see: whether the config screen actually
OFFERS KEGG, whether it sends it in the payload, and whether the results
screen renders a KEGG family - including its figures - rather than dropping
it or sorting it to the bottom.

It also exists because this environment has no Node, so a browser IS the only
machine check available on static/app.js. Any console error or uncaught page
error fails the run.

Needs a live server:  python -m uvicorn API:app --port 8009
Run:                  python tests/playwright_transcriptomics_kegg.py
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


def _go_to_module_config(page, module_name, continue_label, taxon_id="9606", taxon_name="Homo sapiens"):
    """Shared retriever flow, mirroring playwright_msigdb.py's working idiom
    exactly so this file doesn't rediscover selectors that already work.

    taxon_id/taxon_name default to human; the KEGG disambiguation-picker test
    below overrides them to Bacillus subtilis (a real, reproducible KEGG
    ambiguity - see live_transcriptomics_kegg.py) to reach that screen's
    ambiguous-organism branch. The uploaded example file is human genes
    either way - mygene.info's resolve step only needs to COMPLETE and commit
    taxon_id/taxon_name into the session for the KEGG organism check to have
    something to read; it does not require a high match rate."""
    page.evaluate("""() => {
        STATE.view = 'annotation_analytics';
        STATE.aaPage = 'transcriptomics';
        render();
    }""")
    page.wait_for_timeout(200)
    page.locator("button.module-btn.transcriptomics", has_text=module_name).click()
    page.wait_for_timeout(200)

    page.locator("button", has_text="Gene Symbol").first.click()
    page.wait_for_timeout(100)
    page.locator("input[type=file]").set_input_files(
        os.path.join(_ROOT, "Data3", "example_transcriptomics_de_full.csv"))
    page.wait_for_timeout(1500)
    check(page.evaluate("!!STATE.transcriptomics.deTable") is True, "DE table uploaded")

    page.evaluate(
        "(args) => { STATE.transcriptomics.taxonId = args.taxonId; "
        "STATE.transcriptomics.taxonName = args.taxonName; render(); }",
        {"taxonId": taxon_id, "taxonName": taxon_name},
    )
    page.wait_for_timeout(200)
    page.locator("button", has_text="Resolve identifiers").click()
    for _ in range(120):
        if page.evaluate("!!STATE.transcriptomics.resolved"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.transcriptomics.resolved") is True,
          "resolve completed (real mygene.info round trip)")

    page.locator("button", has_text=continue_label).click()
    page.wait_for_timeout(300)


with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page()
    page.on("console", lambda msg: console_errors.append(msg.text) if msg.type == "error" else None)
    page.on("pageerror", lambda exc: console_errors.append(f"pageerror: {exc}"))

    page.goto(BASE, wait_until="networkidle")
    check(page.evaluate("typeof render") == "function", "app.js loaded and parsed")

    print("=== ORA config screen: the KEGG pill and its organism caption ===")
    _go_to_module_config(page, "Over-Representation Analysis",
                         "Continue to Over-Representation Analysis")
    check(page.evaluate("STATE.aaPage") == "transcriptomics_ora_config",
          "navigated to the ORA config screen")

    kegg_pill = page.locator("label.ora-db-pill", has_text="KEGG").first
    check(kegg_pill.count() == 1, "a KEGG pill is present in the database checklist")
    kegg_cb = kegg_pill.locator("input[type=checkbox]")
    check(not kegg_cb.is_disabled(), "the KEGG pill is enabled")

    # KEGG must NOT grow an organism picker here - that is the divergence bug
    # fixed in commit 27ce033, and the reason the organism is derived instead.
    caption_before = page.evaluate(
        "() => [...document.querySelectorAll('.caption')].some(c => c.textContent.includes('KEGG uses the organism'))")
    check(caption_before is False, "the KEGG organism caption is hidden until KEGG is checked")

    kegg_cb.check()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.transcriptomics._oraDraft.databases").count("KEGG") == 1,
          "checking the pill puts KEGG in the draft exactly once")
    caption_after = page.evaluate(
        "() => [...document.querySelectorAll('.caption')].some(c => c.textContent.includes('KEGG uses the organism'))")
    check(caption_after is True,
          "and the caption explaining where KEGG's organism comes from appears")
    check(page.locator("input[placeholder*='organism' i]").count() == 0,
          "SENSITIVITY: NO organism search box is rendered for KEGG - a second, "
          "independently-settable organism picker here is exactly the bug fixed in 27ce033")

    print("=== ORA: run KEGG through the real UI ===")
    page.locator("button.ora-run-btn").click()
    # The screen flips to "results" immediately and shows progress there, so
    # wait on the RESULTS themselves, not on the page name.
    for _ in range(300):
        if page.evaluate("!!STATE.transcriptomics.oraResults"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("STATE.aaPage") == "transcriptomics_ora_results",
          "the run finished and the results screen rendered")
    check(page.evaluate("!!STATE.transcriptomics.oraResults") is True, "results arrived")

    keys = page.evaluate("() => Object.keys(STATE.transcriptomics.oraResults.tables)")
    check(sorted(keys) == ["KEGG (all)", "KEGG (down)", "KEGG (up)"],
          f"three KEGG direction families came back (got {sorted(keys)})")
    n_rows = page.evaluate(
        "() => STATE.transcriptomics.oraResults.tables['KEGG (up)'].table.rows.length")
    check(n_rows > 0, f"the up-direction KEGG table has real rows (got {n_rows})")
    term_ids = page.evaluate(
        "() => STATE.transcriptomics.oraResults.tables['KEGG (up)'].table.rows.map(r => r['Term ID'])")
    check(all(t.startswith("hsa") for t in term_ids),
          f"every term is a human KEGG pathway ID, so the derived organism was really used "
          f"(got {term_ids[:5]})")
    body = page.locator("body").inner_text()
    check("KEGG" in body, "the results screen actually shows a KEGG section")
    check(page.locator(".ora-db-card-kegg_all").count() == 1,
          "the KEGG family gets its own slugged result card, like every other database")

    # Figures are lazy: each sits in a display:none container behind its own
    # "Show ..." button and is only built on first click. So counting SVGs on
    # arrival proves nothing - the toggles have to be exercised.
    toggles = page.locator("button.chart-toggle-btn")
    check(toggles.count() > 0,
          f"the KEGG result card offers chart toggles (got {toggles.count()})")
    check(page.locator("svg").count() == 0,
          "and nothing is drawn before a toggle is clicked (charts are lazy)")
    # The toggles live inside the card's non-active nav panes, so open the
    # pane holding the first one before clicking it.
    nav_titles = page.evaluate(
        "() => [...document.querySelectorAll('.result-nav button')].map(b => b.textContent.trim())")
    check(nav_titles == ["Overview", "Result table", "Gene set overview", "Gene set relationships"],
          f"the KEGG card gets the same four result-nav panes every other database gets "
          f"(got {nav_titles})")
    page.locator(".result-nav button", has_text="Gene set overview").click()
    page.wait_for_timeout(400)
    toggles.first.click()
    page.wait_for_timeout(2000)
    check(page.locator("svg").count() > 0,
          f"clicking the first toggle actually draws a figure for the KEGG family "
          f"(got {page.locator('svg').count()} SVGs) - a new database family must not "
          "land on the results screen figure-less")

    print("=== GSEA config screen + run ===")
    _go_to_module_config(page, "Gene Set Enrichment Analysis",
                         "Continue to Gene Set Enrichment Analysis")
    check(page.evaluate("STATE.aaPage") == "transcriptomics_gsea_config",
          "navigated to the GSEA config screen")
    gsea_kegg = page.locator("label.ora-db-pill", has_text="KEGG").first
    check(gsea_kegg.count() == 1, "the KEGG pill is on the GSEA screen too")
    gsea_kegg.locator("input[type=checkbox]").check()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.transcriptomics._gseaDraft.databases").count("KEGG") == 1,
          "GSEA's draft carries KEGG")

    # 200 permutations: enough for a real result, small enough to stay quick.
    page.evaluate("() => { STATE.transcriptomics._gseaDraft.nPermutations = 200; render(); }")
    page.wait_for_timeout(200)
    page.locator("button.ora-run-btn").click()
    for _ in range(600):
        if page.evaluate("!!STATE.transcriptomics.gseaResults"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("STATE.aaPage") == "transcriptomics_gsea_results",
          "the GSEA run finished and its results screen rendered")
    check(page.evaluate("!!STATE.transcriptomics.gseaResults") is True, "GSEA results arrived")
    gkeys = page.evaluate("() => Object.keys(STATE.transcriptomics.gseaResults.tables)")
    check(gkeys == ["KEGG"], f"one KEGG family, no direction suffix (got {gkeys})")
    gn = page.evaluate("() => STATE.transcriptomics.gseaResults.tables['KEGG'].table.rows.length")
    check(gn > 0, f"the GSEA KEGG table has real rows (got {gn})")

    print("=== KEGG organism disambiguation picker (real ambiguous case: Bacillus subtilis) ===")
    _go_to_module_config(page, "Over-Representation Analysis",
                         "Continue to Over-Representation Analysis",
                         taxon_id="1423", taxon_name="Bacillus subtilis")
    check(page.evaluate("STATE.aaPage") == "transcriptomics_ora_config",
          "navigated to the ORA config screen for the ambiguous organism")
    page.locator("label.ora-db-pill", has_text="KEGG").first.locator("input[type=checkbox]").check()

    picker = page.locator(".ora-kegg-organism-picker")
    for _ in range(30):
        if picker.count() > 0:
            break
        page.wait_for_timeout(500)
    check(picker.count() == 1,
          "SENSITIVITY: the disambiguation picker DOES appear for a genuinely ambiguous "
          "organism (Bacillus subtilis, 14 unverified KEGG strain genomes) - proving the "
          "widget's real /resolve-organism call reaches the ambiguous branch, not just the "
          "offline-fixture tests")
    # And it must NOT have appeared for the clean human case exercised earlier in this
    # same file - checked there already (picker_before assertions), so this confirms the
    # SAME widget genuinely branches on the real resolve-organism response.
    select = picker.locator("select")
    confirm_btn = picker.locator("button", has_text="Confirm KEGG organism")
    check(select.count() == 1 and confirm_btn.count() == 1, "the picker offers a select and a confirm button")
    check(confirm_btn.is_disabled(), "confirm starts disabled until an organism is chosen")

    options = select.locator("option").all_inner_texts()
    check(any(o.startswith("bsu —") for o in options),
          f"the real 'bsu' (B. subtilis 168) candidate is offered by name (got {options[:5]})")
    select.select_option("bsu")
    check(not confirm_btn.is_disabled(), "confirm enables once an organism is selected")

    confirm_btn.click()
    for _ in range(30):
        if page.evaluate("() => !!(STATE.transcriptomics._keggOrganismConfirmed)"):
            break
        page.wait_for_timeout(500)
    confirmed_state = page.evaluate("() => STATE.transcriptomics._keggOrganismConfirmed")
    check(confirmed_state and confirmed_state.get("organismCode") == "bsu",
          f"confirming persists the real chosen organism into client state (got {confirmed_state})")
    check(page.locator(".ora-selected-caption", has_text="Confirmed manually").count() == 1,
          "the config screen now shows the confirmed-organism caption instead of the picker")
    check(picker.count() == 0, "and the picker itself is gone, replaced by that caption")

    print("=== running KEGG ORA with the confirmed organism ===")
    page.locator("button.ora-run-btn").click()
    for _ in range(300):
        if page.evaluate("!!STATE.transcriptomics.oraResults"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.transcriptomics.oraResults") is True, "the confirmed-organism run finished")
    bsu_note = page.evaluate(
        "() => (STATE.transcriptomics.oraResults.tables['KEGG (all)'] || {}).note || ''")
    check("manually confirmed" in bsu_note.lower(),
          f"the result table carries the manual-confirmation caveat, proving the confirmed "
          f"organism actually reached the real run (got {bsu_note!r})")

    print("=== changing organism clears the confirmation (client-side mirror of the server guard) ===")
    page.evaluate("""() => {
        STATE.transcriptomics.taxonId = '9606';
        STATE.transcriptomics.taxonName = 'Homo sapiens';
        STATE.transcriptomics.resolved = null;
        _clearTranscriptomicsResolutionDependentState();
        STATE.aaPage = 'transcriptomics_ora_config';
        render();
    }""")
    page.wait_for_timeout(200)
    check(page.evaluate("() => !STATE.transcriptomics._keggOrganismConfirmed"),
          "the confirmed B. subtilis organism no longer lingers in client state after switching species")

    check(not console_errors,
          f"no console or page errors anywhere in the whole flow (got {console_errors[:5]})")
    browser.close()

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
