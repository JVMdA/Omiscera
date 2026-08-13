"""Regression test for a real, user-reported bug in buildTable() (static/
app.js), the SHARED table renderer used by every results screen in the app
(Genomics ORA/Burden/Annotation, Transcriptomics ORA/GSEA/Regulon) - not
MSigDB-specific code, even though MSigDB's large Hallmark gene sets are what
first exposed it.

THE BUG. A comma-joined gene-list cell (Overlapping Genes / Leading Edge
Genes / Genes) can hold hundreds of symbols. Wrapped inside a narrow column,
that ONE cell's text can inflate its WHOLE ROW to thousands of pixels tall -
confirmed live: an MSigDB Hallmark term with 185 overlapping genes produced a
row 3295px tall. Table cells default to vertical-align:middle, so every
OTHER (short) cell in that same row - Term ID, k, K, FDR, everything a user
actually scans - gets vertically centered ~1600px down, far outside the
results table's 480px scroll viewport. The table LOOKS completely empty
(just a solid colored strip, no visible text) even though the data is
correctly in the DOM.

THE FIX has two independent layers, both tested here: (1) long-list columns
are capped to 3 lines via CSS line-clamp, so they can no longer cause the row
explosion in the first place; (2) vertical-align:top on every cell, as a
second, independent safeguard for any OTHER column that turns out to hold
unexpectedly long content in the future.

Needs a live server (to load app.js and its CSS) - no MSigDB network calls,
purely synthetic fixtures.
"""
import sys

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

    print("=== a long 'Overlapping Genes' cell no longer explodes its row's height ===")
    long_gene_list = ", ".join(f"GENE{i}" for i in range(200))  # ~1400 chars, similar order to the real repro
    result = page.evaluate(
        """(geneList) => {
        const div = document.createElement('div');
        document.body.appendChild(div);
        const payload = {
            columns: ['Term ID', 'Term Name', 'k', 'K', 'FDR', 'Overlapping Genes'],
            rows: [
                { 'Term ID': 'REAL_TERM', 'Term Name': 'A term with a huge overlap', 'k': 185, 'K': 188, 'FDR': 1e-280,
                  'Overlapping Genes': geneList },
            ],
        };
        const wrap = buildTable(payload, { maxHeight: 480 });
        div.appendChild(wrap);
        const tr = div.querySelector('tbody tr');
        const cells = [...tr.querySelectorAll('td')];
        return {
            rowHeight: tr.getBoundingClientRect().height,
            firstCellText: cells[0].textContent,
            firstCellTop: cells[0].getBoundingClientRect().top,
            trTop: tr.getBoundingClientRect().top,
            longCellHasClass: cells[5].classList.contains('td-long-list'),
            longCellTitle: cells[5].getAttribute('title'),
            longCellVisibleText: cells[5].querySelector('.td-long-list-content').textContent,
        };
    }""",
        long_gene_list,
    )
    check(result["rowHeight"] < 200,
          f"row height stays reasonable (<200px) despite a ~1400-char gene list "
          f"(got {result['rowHeight']:.1f}px - the real bug produced 3295px)")
    check(abs(result["firstCellTop"] - result["trTop"]) < 5,
          f"the FIRST cell ('Term ID') renders at the TOP of its row, not vertically centered deep inside it "
          f"(cell top={result['firstCellTop']:.1f}, row top={result['trTop']:.1f})")
    check(result["firstCellText"] == "REAL_TERM",
          f"the short cell's text is genuinely present and correct (got {result['firstCellText']!r})")
    check(result["longCellHasClass"] is True, "the long-list cell gets the capping CSS class")
    check(result["longCellTitle"] == long_gene_list,
          "the FULL, untruncated gene list is preserved in the title tooltip - this is a DISPLAY cap only")
    check(result["longCellVisibleText"] == long_gene_list,
          "the full text is still in the DOM (line-clamp is visual-only, not a data truncation) - "
          "a user can still Ctrl+F / select-all it even though only 3 lines are visible")

    print("=== a SHORT gene list is rendered plainly, no unnecessary wrapper ===")
    short_result = page.evaluate("""() => {
        const div = document.createElement('div');
        document.body.appendChild(div);
        const payload = {
            columns: ['Term ID', 'Overlapping Genes'],
            rows: [{ 'Term ID': 'SHORT_TERM', 'Overlapping Genes': 'GENE1, GENE2, GENE3' }],
        };
        const wrap = buildTable(payload, {});
        div.appendChild(wrap);
        const cell = div.querySelectorAll('tbody tr td')[1];
        return { hasClass: cell.classList.contains('td-long-list'), text: cell.textContent };
    }""")
    check(short_result["hasClass"] is False,
          f"a short (<=100 char) gene list is NOT wrapped in the capping class - no unnecessary DOM nesting "
          f"for the common case (got hasClass={short_result['hasClass']})")
    check(short_result["text"] == "GENE1, GENE2, GENE3", "and its text renders normally")

    print("=== SENSITIVITY: the row-height check is not vacuous ===")
    # Re-run the exact same measurement with the fix's CSS class deliberately
    # withheld (simulating the pre-fix buildTable, which never special-cased
    # long columns at all) and confirm the row genuinely DOES explode.
    broken_result = page.evaluate(
        """(geneList) => {
        const div = document.createElement('div');
        document.body.appendChild(div);
        // Mirrors the PRE-FIX buildTable: every cell is a plain <td>, no
        // long-list capping, relying on the browser's default wrapping.
        const table = document.createElement('table');
        table.className = 'data-table';
        const tbody = document.createElement('tbody');
        const tr = document.createElement('tr');
        ['REAL_TERM', geneList].forEach((v) => {
            const td = document.createElement('td');
            td.textContent = v;
            tr.appendChild(td);
        });
        tbody.appendChild(tr);
        table.appendChild(tbody);
        const wrap = document.createElement('div');
        wrap.className = 'table-wrap';
        const scroller = document.createElement('div');
        scroller.className = 'table-scroll';
        scroller.style.maxHeight = '480px';
        scroller.appendChild(table);
        wrap.appendChild(scroller);
        div.appendChild(wrap);
        return { rowHeight: tr.getBoundingClientRect().height };
    }""",
        long_gene_list,
    )
    # Relative, not an absolute magic number: the exact pixel height depends
    # on the test browser's viewport width (a wider window wraps the same
    # text across fewer, longer lines) - but the FIXED cell must always be
    # dramatically shorter than the unfixed one for the same input, in any
    # width. The live repro measured 3295px broken vs ~75px fixed (~44x);
    # requiring merely 2x here is a conservative floor that still fails
    # if someone quietly regresses the cap back toward the old behaviour.
    check(broken_result["rowHeight"] > result["rowHeight"] * 2,
          f"the PRE-FIX rendering (no capping) produces a row at least 2x taller than the fixed version for "
          f"the exact same input (broken={broken_result['rowHeight']:.1f}px vs fixed={result['rowHeight']:.1f}px) "
          f"- proving the check above is not vacuous")

    check(len(errors) == 0, f"no page errors while exercising buildTable (got {errors})")
    browser.close()

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
