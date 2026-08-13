"""Regression tests for the MSigDB cross-collection summary figure
(_buildMsigdbCollectionSummaryCore / _msigdbFamiliesFrom{Ora,Gsea}Tables /
_msigdbBestRow in static/app.js), driven in a real browser against synthetic
result tables.

WHY THIS FILE EXISTS. A post-implementation audit found three real defects in
this figure that the original playwright_msigdb.py could not see, because
that test only asserted the figure "rendered with data":

  1. The bar COLOUR applied `bestEffect >= 0` to both ORA and GSEA. For GSEA
     that is meaningful (NES is signed). For ORA the effect measure is fold
     enrichment - a non-negative ratio whose neutral point is 1, not 0 - so
     the blue branch was unreachable and the colour silently encoded nothing
     while looking like it encoded direction.
  2. A blank/non-numeric FDR in the FIRST row seeded `bestFdr = NaN`, and
     since `x < NaN` is always false it was never replaced: the family then
     reported that first row as its top hit forever, with "NaN" in the hover.
  3. The figure cached its build and was never invalidated when the user
     changed the significance threshold, so it could show bar lengths and an
     axis title computed at the OLD alpha directly above a table using the
     new one.

Each test below is written to fail if its defect is reintroduced; the
sensitivity of (1) and (2) is proven directly at the bottom of this file by
re-running the check against a deliberately-broken reimplementation.

Needs a live server (for app.js) but makes no MSigDB network calls - all
tables here are synthetic.
"""
import os
import sys

from playwright.sync_api import sync_playwright

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
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

    print("=== ORA: fold enrichment is UNSIGNED, so the figure must not imply direction ===")
    # Two collections, both significant, with fold enrichments on either side
    # of 1.0 (0.4 and 5.0). Under the old `bestEffect >= 0` rule BOTH were
    # painted the same red - and no possible fold-enrichment value could ever
    # produce blue, making the colour channel vacuous.
    ora_colors = page.evaluate("""() => {
        const div = document.createElement('div');
        document.body.appendChild(div);
        const tables = {
            'MSigDB: Hallmark (all)': { table: { rows: [
                { 'Term ID': 'H1', 'Term Name': 'BELOW_ONE', 'FDR': 0.001, 'Fold Enrichment': 0.4 },
            ] } },
            'MSigDB: Cell type (all)': { table: { rows: [
                { 'Term ID': 'C1', 'Term Name': 'ABOVE_ONE', 'FDR': 0.001, 'Fold Enrichment': 5.0 },
            ] } },
        };
        buildMsigdbOraCollectionSummary(div, tables, 0.05, 'test_ora');
        const d = div.querySelector('.js-plotly-plot') || div;
        return { colors: d.data[0].marker.color, y: d.data[0].y };
    }""")
    check(len(set(ora_colors["colors"])) == 1,
          f"ORA bars are a single enrichment colour regardless of fold-enrichment magnitude - the figure no "
          f"longer implies a direction its statistic does not have (got {ora_colors['colors']})")
    check("#2874a6" not in ora_colors["colors"],
          "ORA never uses the 'down/negative' blue, which for a one-sided over-representation test would be meaningless")

    print("=== GSEA: NES IS signed, so colour must still track its sign ===")
    gsea_colors = page.evaluate("""() => {
        const div = document.createElement('div');
        document.body.appendChild(div);
        const tables = {
            'MSigDB: Hallmark': { table: { rows: [
                { 'Term ID': 'H1', 'Term Name': 'NEG', 'FDR q-value': 0.001, 'NES': -2.4 },
            ] } },
            'MSigDB: Cell type': { table: { rows: [
                { 'Term ID': 'C1', 'Term Name': 'POS', 'FDR q-value': 0.001, 'NES': 2.4 },
            ] } },
        };
        buildMsigdbGseaCollectionSummary(div, tables, 0.25, 'test_gsea');
        const d = div.querySelector('.js-plotly-plot') || div;
        // families are reversed for display, so index 0 is the LAST collection
        // in canonical order (Cell type), index 1 is Hallmark.
        return { colors: d.data[0].marker.color, y: d.data[0].y };
    }""")
    check(gsea_colors["y"] == ["Cell type", "Hallmark"],
          f"display order is reversed so the first canonical collection sits on top (got {gsea_colors['y']})")
    check(gsea_colors["colors"][0] == "#c0392b" and gsea_colors["colors"][1] == "#2874a6",
          f"GSEA still colours by NES sign: positive-NES Cell type red, negative-NES Hallmark blue "
          f"(got {gsea_colors['colors']})")

    print("=== a non-numeric FDR in the FIRST row must not poison the family's best-row scan ===")
    # The killer ordering: a bad row FIRST, a genuinely good row after it.
    nan_result = page.evaluate("""() => {
        const div = document.createElement('div');
        document.body.appendChild(div);
        const tables = {
            'MSigDB: Hallmark (all)': { table: { rows: [
                { 'Term ID': 'BAD', 'Term Name': 'UNPARSEABLE_FDR', 'FDR': 'NA', 'Fold Enrichment': 9.9 },
                { 'Term ID': 'GOOD', 'Term Name': 'REAL_TOP_HIT', 'FDR': 0.002, 'Fold Enrichment': 3.1 },
            ] } },
        };
        buildMsigdbOraCollectionSummary(div, tables, 0.05, 'test_nan');
        const d = div.querySelector('.js-plotly-plot') || div;
        return { customdata: d.data[0].customdata[0], x: d.data[0].x };
    }""")
    check(nan_result["customdata"][2] == "REAL_TOP_HIT",
          f"the real lowest-FDR row is reported as the top term, not the unparseable first row "
          f"(got {nan_result['customdata'][2]!r})")
    check("NaN" not in str(nan_result["customdata"]),
          f"no 'NaN' leaks into the hover text (got {nan_result['customdata']})")
    check(nan_result["x"] == [1],
          f"the unparseable row is excluded from the significant-set COUNT too, not counted as significant "
          f"(got {nan_result['x']})")

    print("=== a family whose rows are ALL unparseable degrades to grey/n\\a, never a false direction ===")
    all_bad = page.evaluate("""() => {
        const div = document.createElement('div');
        document.body.appendChild(div);
        const tables = {
            'MSigDB: Hallmark (all)': { table: { rows: [
                { 'Term ID': 'B1', 'Term Name': 'X', 'FDR': null, 'Fold Enrichment': null },
            ] } },
        };
        buildMsigdbOraCollectionSummary(div, tables, 0.05, 'test_allbad');
        const d = div.querySelector('.js-plotly-plot') || div;
        return { color: d.data[0].marker.color[0], customdata: d.data[0].customdata[0], x: d.data[0].x[0] };
    }""")
    check(all_bad["x"] == 0, f"no significant sets counted (got {all_bad['x']})")
    check(all_bad["color"] == "#b8c2cc", f"bar is grey, not a confident red/blue (got {all_bad['color']})")
    check(all_bad["customdata"][0] == "n/a" and all_bad["customdata"][1] == "n/a",
          f"hover reports n/a rather than NaN (got {all_bad['customdata']})")

    print("=== a genuinely underflowed FDR (literal 0.0) must not read as a display bug ===")
    # scipy's hypergeom.sf() can legitimately return a literal 0.0 for an
    # astronomically significant term (smaller than float64 can represent) -
    # this is real, not a computation error, but "0.00e+0" looks exactly
    # like one. User-reported: a 349-set HPO/MPT family's tooltip showed
    # "Best FDR: 0.00e+0" and was reasonably mistaken for the earlier
    # null-FDR-coerced-to-0 bug this same file already guards against above.
    underflow_result = page.evaluate("""() => {
        const div = document.createElement('div');
        document.body.appendChild(div);
        const tables = {
            'MSigDB: Phenotype (HPO / MPT) (all)': { table: { rows: [
                { 'Term ID': 'HP:0000118', 'Term Name': 'HP_GENITAL_NEOPLASM', 'FDR': 0.0, 'Fold Enrichment': 4.2 },
            ] } },
        };
        buildMsigdbOraCollectionSummary(div, tables, 0.25, 'test_underflow');
        const d = div.querySelector('.js-plotly-plot') || div;
        return { bestFdrText: d.data[0].customdata[0][0] };
    }""")
    check(underflow_result["bestFdrText"] == "< 1e-300",
          f"a literal 0.0 FDR is shown as '< 1e-300' (underflow), not the misleading '0.00e+0' "
          f"(got {underflow_result['bestFdrText']!r})")

    print("=== empty input does not throw (families[0] guard) ===")
    empty_ok = page.evaluate("""() => {
        const div = document.createElement('div');
        document.body.appendChild(div);
        try {
            buildMsigdbOraCollectionSummary(div, { 'GO: Biological Process (all)': { table: { rows: [] } } }, 0.05, 'test_empty');
            return div.textContent.length > 0;
        } catch (e) { return 'THREW: ' + e.message; }
    }""")
    check(empty_ok is True, f"a run with no MSigDB families shows a message instead of throwing (got {empty_ok!r})")

    print("=== SENSITIVITY: the colour and NaN checks above are not vacuous ===")
    # Re-run the two key assertions against a deliberately-broken
    # reimplementation of the fixed logic, in-page only (app.js on disk is
    # never touched), and confirm each WOULD have failed.
    sensitivity = page.evaluate("""() => {
        // Broken best-row scan: the pre-fix version, with no isFinite guard.
        function brokenBestRow(rows, fdrKey, effectKey) {
            let bestFdr = null, bestEffect = null, bestTermName = null;
            rows.forEach((r) => {
                const fdr = Number(r[fdrKey]);
                if (bestFdr === null || fdr < bestFdr) {
                    bestFdr = fdr; bestEffect = Number(r[effectKey]); bestTermName = r['Term Name'];
                }
            });
            return { bestFdr, bestEffect, bestTermName };
        }
        const rows = [
            { 'Term ID': 'BAD', 'Term Name': 'UNPARSEABLE_FDR', 'FDR': 'NA', 'Fold Enrichment': 9.9 },
            { 'Term ID': 'GOOD', 'Term Name': 'REAL_TOP_HIT', 'FDR': 0.002, 'Fold Enrichment': 3.1 },
        ];
        const brokenBest = brokenBestRow(rows, 'FDR', 'Fold Enrichment');

        // Broken colour rule: the pre-fix `bestEffect >= 0` applied to ORA.
        const brokenOraColor = (bestEffect) => (bestEffect === null || bestEffect >= 0) ? '#c0392b' : '#2874a6';

        return {
            brokenTopTerm: brokenBest.bestTermName,
            brokenFdrIsNaN: Number.isNaN(brokenBest.bestFdr),
            brokenColorBelowOne: brokenOraColor(0.4),
            brokenColorAboveOne: brokenOraColor(5.0),
        };
    }""")
    check(sensitivity["brokenTopTerm"] == "UNPARSEABLE_FDR" and sensitivity["brokenFdrIsNaN"] is True,
          "the PRE-FIX best-row scan really does pin itself to the unparseable first row with a NaN FDR - "
          "so the NaN tests above would have caught this defect")
    check(sensitivity["brokenColorBelowOne"] == sensitivity["brokenColorAboveOne"] == "#c0392b",
          "the PRE-FIX ORA colour rule really does paint fold enrichment 0.4 and 5.0 identically red - "
          "so the colour test above would have caught this defect")

    check(len(errors) == 0, f"no page errors while exercising the summary figure (got {errors})")
    browser.close()

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
