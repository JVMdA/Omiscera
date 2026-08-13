"""Regression test for a real, user-reported bug: the Genomics 'Network
Propagation / Module Detection' and 'Gene Function Annotation' results
screens used to append their makeResultNav host straight to the page with
no wrapping card, unlike every other results screen (Transcriptomics
GSEA/ORA, Genomics ORA/Burden), which all sit inside the standard bordered
'.ora-db-card' (blue top accent). Fixed by wrapping navHost in that same
card class - see the identical fix (and its test) for the Transcriptomics
TF/Regulon screen in playwright_regulon.py.

This is also, as far as this suite goes, the FIRST browser coverage either
of these two Genomics screens has ever had - there was no genomics
Playwright file before this one. Drives the REAL upload -> resolve ->
configure -> run flow for each (no synthetic STATE injection - a hand-built
fixture risks silently drifting from the real payload shape and giving
false confidence).

Needs a live server and real network access (UniProt resolve, STRING
network fetch, GO/UniProt annotation lookups) - not mocked, takes a
minute or two.
"""
import os
import sys

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


def card_check(page):
    return page.evaluate("""() => {
        const nav = document.querySelector('.result-nav-layout');
        return !!(nav && nav.closest('.ora-db-card'));
    }""")


with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page()
    errors = []
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    page.goto(BASE, wait_until="networkidle")
    page.wait_for_timeout(300)

    print("=== Network Propagation / Module Detection ===")
    page.evaluate("""() => {
        STATE.view = 'annotation_analytics';
        STATE.aaPage = 'genomics';
        render();
    }""")
    page.wait_for_timeout(200)
    page.locator("button.module-btn.genomics", has_text="Network Propagation / Module Detection").click()
    page.wait_for_timeout(200)
    page.locator("button", has_text="Gene Symbol").first.click()
    page.wait_for_timeout(150)

    # Organism set directly - this test is about the results-screen card,
    # not the STRING species-search widget, which has no bearing on that.
    page.evaluate("""() => { STATE.ora.taxonId = '9606'; STATE.ora.taxonName = 'Homo sapiens'; render(); }""")
    page.wait_for_timeout(150)

    page.locator("input[type=file]").set_input_files(
        os.path.join(_ROOT, "Data2", "example_netprop_input_symbol.csv")
    )
    page.wait_for_timeout(1000)

    resolve_btn = page.locator("button", has_text="Resolve genes")
    check(resolve_btn.count() >= 1 and not resolve_btn.first.is_disabled(), "netprop: Resolve genes button enabled")
    resolve_btn.first.click()
    for _ in range(120):
        if page.evaluate("!!STATE.ora.resolved"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.ora.resolved") is True, "netprop: resolve completed (real UniProt round trip)")
    page.wait_for_timeout(200)
    if page.evaluate("!STATE.ora.finalTable && (STATE.ora.resolved.ambiguous||[]).length > 0"):
        page.locator("button", has_text="Discard all unresolved").click()
        page.wait_for_timeout(300)

    continue_btn = page.locator("button", has_text="Continue to analysis")
    check(continue_btn.count() >= 1, "netprop: 'Continue to analysis' button present")
    continue_btn.first.click()
    page.wait_for_timeout(300)

    run_btn = page.locator("button.ora-run-btn")
    check(run_btn.count() >= 1 and not run_btn.first.is_disabled(), "netprop: Run button enabled")
    run_btn.first.click()
    for _ in range(180):
        if page.evaluate("!!STATE.ora.netpropResults"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.ora.netpropResults") is True, "netprop: job completed")
    page.wait_for_timeout(400)

    check(card_check(page), "netprop: results nav sits inside '.ora-db-card' (was bare on the page before the fix)")

    print("=== Gene Function Annotation ===")
    page.evaluate("""() => {
        STATE.view = 'annotation_analytics';
        STATE.aaPage = 'genomics';
        render();
    }""")
    page.wait_for_timeout(200)
    page.locator("button.module-btn.genomics", has_text="Gene Function Annotation").click()
    page.wait_for_timeout(200)
    page.locator("button", has_text="Gene Symbol").first.click()
    page.wait_for_timeout(150)
    page.evaluate("""() => { STATE.ora.taxonId = '9606'; STATE.ora.taxonName = 'Homo sapiens'; render(); }""")
    page.wait_for_timeout(150)

    # Small, fast hand-built list - annotation is descriptive (no stats), so
    # a tiny gene set is fine and keeps this quick.
    page.evaluate("""() => {
        const dt = new DataTransfer();
        const csv = 'Gene Symbol\\nTP53\\nBRCA1\\nMYC\\n';
        const file = new File([csv], 'tiny_genes.csv', { type: 'text/csv' });
        dt.items.add(file);
        const input = document.querySelector('input[type=file]');
        input.files = dt.files;
        input.dispatchEvent(new Event('change', { bubbles: true }));
    }""")
    page.wait_for_timeout(1200)

    resolve_btn = page.locator("button", has_text="Resolve genes")
    check(resolve_btn.count() >= 1 and not resolve_btn.first.is_disabled(), "annotation: Resolve genes button enabled")
    resolve_btn.first.click()
    for _ in range(120):
        if page.evaluate("!!STATE.ora.resolved"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.ora.resolved") is True, "annotation: resolve completed (real UniProt round trip)")
    page.wait_for_timeout(200)
    if page.evaluate("!STATE.ora.finalTable && (STATE.ora.resolved.ambiguous||[]).length > 0"):
        page.locator("button", has_text="Discard all unresolved").click()
        page.wait_for_timeout(300)

    continue_btn = page.locator("button", has_text="Continue to analysis")
    check(continue_btn.count() >= 1, "annotation: 'Continue to analysis' button present")
    continue_btn.first.click()
    page.wait_for_timeout(300)

    run_btn = page.locator("button.ora-run-btn")
    check(run_btn.count() >= 1 and not run_btn.first.is_disabled(), "annotation: Run button enabled")
    run_btn.first.click()
    for _ in range(120):
        if page.evaluate("!!STATE.ora.annotationResults"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.ora.annotationResults") is True, "annotation: job completed")
    page.wait_for_timeout(400)

    check(card_check(page), "annotation: results nav sits inside '.ora-db-card' (was bare on the page before the fix)")

    check(len(errors) == 0, f"no page errors during either flow (got {errors})")
    browser.close()

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
