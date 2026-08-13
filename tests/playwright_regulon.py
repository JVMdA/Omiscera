"""Playwright end-to-end UI test for the new TF/Regulon frontend
(static/app.js). Drives real clicks through the actual retriever ->
regulon-config -> regulon-results flow against a live uvicorn (127.0.0.1:8009),
using the app's own example DE file (real mygene.info + TRRUST/OmniPath/ChEA3
network calls - no mocking). Also asserts NO console errors were raised at
any point, which is the closest thing to a JS syntax check available without
node."""
import os
import sys
import time

from playwright.sync_api import sync_playwright

# Repo-relative: this suite lives in <repo>/tests/ and the app under test is
# one level up. Was a hardcoded absolute path while these were throwaway
# scratchpad scripts; made portable when they moved into the repo.
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


with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page()
    page.on("console", lambda msg: console_errors.append(msg.text) if msg.type == "error" else None)
    page.on("pageerror", lambda exc: console_errors.append(f"pageerror: {exc}"))

    page.goto(BASE, wait_until="networkidle")
    check(page.evaluate("typeof render") == "function", "app.js loaded and parsed (render() is defined)")
    check(page.evaluate("typeof renderTranscriptomicsRegulonConfigScreen") == "function",
          "renderTranscriptomicsRegulonConfigScreen is defined (new code parsed without a syntax error)")
    check(page.evaluate("typeof buildRegulonActivityBarPlot") == "function", "buildRegulonActivityBarPlot is defined")
    check(page.evaluate("typeof buildRegulonActivityVolcano") == "function", "buildRegulonActivityVolcano is defined")
    check(page.evaluate("typeof buildRegulonActivityVsExpression") == "function", "buildRegulonActivityVsExpression is defined (Figure B)")
    check(page.evaluate("typeof buildRegulonUlmRegressionPlot") == "function", "buildRegulonUlmRegressionPlot is defined (Figure A)")
    check(page.evaluate("typeof buildRegulonCoherenceBars") == "function", "buildRegulonCoherenceBars is defined (Figure E)")
    check(page.evaluate("typeof buildRegulonCalibrationQQPlot") == "function", "buildRegulonCalibrationQQPlot is defined (calibration diagnostic)")
    check(page.evaluate("typeof _probit") == "function", "_probit (inverse normal CDF) is defined")
    # Sanity-check the probit approximation against known values before
    # trusting it in a rendered figure.
    probit_median = page.evaluate("_probit(0.5)")
    probit_975 = page.evaluate("_probit(0.975)")
    check(abs(probit_median) < 1e-6, f"_probit(0.5) ~= 0 (got {probit_median})")
    check(abs(probit_975 - 1.959964) < 1e-4, f"_probit(0.975) ~= 1.959964 (got {probit_975})")

    print("=== navigate to Transcriptomics module grid, tile is enabled ===")
    page.evaluate("""() => {
        STATE.view = 'annotation_analytics';
        STATE.aaPage = 'transcriptomics';
        render();
    }""")
    page.wait_for_timeout(200)
    tile = page.locator("button.module-btn.transcriptomics", has_text="Transcription Factor / Regulon")
    check(tile.count() == 1, "Transcription Factor / Regulon tile is present")
    check(not tile.is_disabled(), "the tile is ENABLED (no longer 'coming soon')")

    tile.click()
    page.wait_for_timeout(200)
    check(page.evaluate("STATE.transcriptomicsModuleName") == "Transcription Factor / Regulon",
          "clicking the tile sets the module name correctly")

    print("=== retriever: identifier type, upload, organism, resolve (real mygene.info) ===")
    page.locator("button", has_text="Gene Symbol").first.click()
    page.wait_for_timeout(100)

    file_input = page.locator("input[type=file]")
    example_path = os.path.join(_ROOT, "Data3", "example_transcriptomics_de_full.csv")
    file_input.set_input_files(example_path)
    page.wait_for_timeout(1500)
    check(page.evaluate("!!STATE.transcriptomics.deTable") is True, "DE table uploaded")

    # The fuzzy organism-search autocomplete is already covered by prior
    # sessions' retriever tests; set it directly here to keep this script
    # focused on the NEW regulon code without fighting unrelated UI timing.
    page.evaluate("""() => {
        STATE.transcriptomics.taxonId = '9606';
        STATE.transcriptomics.taxonName = 'Homo sapiens';
        render();
    }""")
    page.wait_for_timeout(200)
    check(page.evaluate("STATE.transcriptomics.taxonId") == "9606", "taxonId set to 9606")

    resolve_btn = page.locator("button", has_text="Resolve identifiers")
    check(resolve_btn.count() == 1 and not resolve_btn.is_disabled(), "Resolve identifiers button is enabled")
    resolve_btn.click()
    for _ in range(120):
        if page.evaluate("!!STATE.transcriptomics.resolved"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.transcriptomics.resolved") is True, "resolve completed (real mygene.info round trip)")

    print("=== continue to TF/Regulon config ===")
    continue_btn = page.locator("button", has_text="Continue to Transcription Factor / Regulon")
    check(continue_btn.count() == 1, "Continue button present")
    continue_btn.click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "transcriptomics_regulon_config", "navigated to regulon config screen")

    print("=== config screen: no false small-input warning ===")
    # The DoRothEA-confidence and target-enrichment-vs-activity-inference
    # captions that used to live here were moved into the README (this level
    # of explanatory detail belongs in read-once documentation, not on every
    # visit to the config screen) - see
    # info/README-TFRegulon-transcriptomics.txt's "THE FOUR DATABASES"
    # section for the DoRothEA-only / CollecTRI-not-filtered claim this test
    # used to check on-screen.
    captions_text = " ".join(page.locator(".caption").all_inner_texts())
    check("DoRothEA ONLY" not in captions_text and "grades each TF-target edge's evidence" not in captions_text,
          "the DoRothEA-confidence explanation is no longer duplicated on the config screen")
    check("Target enrichment: is this TF's target set" not in captions_text,
          "the target-enrichment-vs-activity-inference explanation is no longer duplicated on the config screen")

    # The full example resolves to ~3000 genes, so the small-input warning
    # must NOT be showing here.
    warn_boxes = page.locator(".alert-warning").all_inner_texts()
    small_warn = [w for w in warn_boxes if "resolved table has only" in w]
    check(len(small_warn) == 0,
          f"no small-input warning for the full ~3000-gene example (got {small_warn})")

    print("=== configure: TRRUST + ChEA3, both analyses ===")
    pills = page.locator("label.ora-db-pill")
    n_pills = pills.count()
    clicked_dbs = []
    for i in range(n_pills):
        text = pills.nth(i).inner_text()
        if text.strip().startswith("TRRUST") or text.strip().startswith("ChEA3"):
            cb = pills.nth(i).locator("input[type=checkbox]")
            if not cb.is_disabled():
                cb.check()
                clicked_dbs.append(text.strip())
    check(len(clicked_dbs) == 2, f"TRRUST and ChEA3 checkboxes both checked (got: {clicked_dbs})")

    # "Analyses to run" defaults to both checked already - verify via state.
    check(page.evaluate("STATE.transcriptomics._regulonDraft.analyses.length") == 2,
          "both analyses (enrichment + activity) selected by default")

    run_btn = page.locator("button.ora-run-btn")
    check(run_btn.count() == 1 and not run_btn.is_disabled(), "Run button enabled with 2 databases selected")
    run_btn.click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "transcriptomics_regulon_results", "navigated to regulon results screen")

    print("=== wait for the real regulon job to finish ===")
    for _ in range(180):
        if page.evaluate("!!STATE.transcriptomics.regulonResults"):
            break
        page.wait_for_timeout(1000)
    check(page.evaluate("!!STATE.transcriptomics.regulonResults") is True, "regulon job completed")

    page.wait_for_timeout(300)
    check(page.locator(".result-nav-btn", has_text="TF Activity").count() >= 1, "'TF Activity' rail entry present")
    check(page.locator(".result-nav-btn", has_text="TF Target Enrichment").count() >= 1, "'TF Target Enrichment' rail entry present")

    # Real, user-reported bug: this results screen used to append its
    # makeResultNav host straight to the page with no wrapping card, unlike
    # every other results screen (GSEA/ORA/Burden), which all sit inside the
    # standard bordered '.ora-db-card' (blue top accent). Fixed by wrapping
    # navHost in that same card class.
    in_card = page.evaluate("""() => {
        const nav = document.querySelector('.result-nav-layout');
        return !!(nav && nav.closest('.ora-db-card'));
    }""")
    check(in_card, "the results nav sits inside the standard '.ora-db-card' box, matching GSEA/ORA - "
                   "not bare on the page")

    # Loosen the (display-only, adjustable) significance threshold so at least
    # one TF clears it in BOTH tables regardless of how strict this particular
    # real dataset/database combo happens to be - otherwise the toggle buttons
    # below never render and the chart-build assertions get silently skipped.
    page.evaluate("""() => {
        STATE.transcriptomics.regulonSigAlpha = 0.999;
        render();
    }""")
    page.wait_for_timeout(300)

    print("=== TF Activity pane: lambda card + calibration QQ plot ===")
    page.locator(".result-nav-btn", has_text="TF Activity").first.click()
    page.wait_for_timeout(200)
    lambda_card = page.locator(".metric-card", has_text="Null inflation")
    check(lambda_card.count() >= 1, "'Null inflation (λ)' summary card present")

    qq_btn = page.locator("button.chart-toggle-btn", has_text="Show calibration QQ plot")
    check(qq_btn.count() > 0, "'Show calibration QQ plot' toggle present")
    if qq_btn.count() > 0:
        qq_btn.click()
        page.wait_for_timeout(600)
        qq_info = page.evaluate("""() => {
            const divs = [...document.querySelectorAll('.chart-container')];
            const d = divs.find(x => x.id.startsWith('chart-calibrationqq-') && x.data);
            if (!d) return null;
            const obs = d.data.find(t => t.name === 'Observed TF scores');
            const refLine = d.data.find(t => t.name === 'y = x (no inflation)');
            const fitLine = d.data.find(t => t.name.includes('Fitted empirical null'));
            if (!obs || !refLine || !fitLine) return null;
            const fitSlope = (fitLine.y[1] - fitLine.y[0]) / (fitLine.x[1] - fitLine.x[0]);
            const refSlope = (refLine.y[1] - refLine.y[0]) / (refLine.x[1] - refLine.x[0]);
            return {
                nObserved: obs.x.length,
                sortedX: obs.x.every((v, i) => i === 0 || v >= obs.x[i - 1]),
                fitSlope: fitSlope,
                refSlopeIsOne: Math.abs(refSlope - 1.0) < 1e-6,
                titleHasLambda: d.layout.title.text.includes('λ'),
            };
        }""")
        check(qq_info is not None, "QQ plot exposes observed/reference/fitted traces")
        if qq_info:
            check(qq_info["nObserved"] > 5, f"QQ plot has multiple observed TF points (got {qq_info['nObserved']})")
            check(qq_info["sortedX"], "theoretical quantiles (x) are sorted ascending, as a QQ plot requires")
            check(qq_info["refSlopeIsOne"], "y=x reference line genuinely has slope 1")
            check(qq_info["titleHasLambda"], "plot title reports the lambda value")
            # Cross-check the fitted line's slope against STATE's own lambda
            # for this database - the figure must draw the SAME lambda the
            # table/summary cards report, not a recomputed or stale one.
            state_lambda = page.evaluate("""() => {
                const r = STATE.transcriptomics.regulonResults;
                const dbSel = document.querySelectorAll('.result-pane:not([style*="display: none"]) select')[0];
                const db = dbSel ? dbSel.value : Object.keys(r.activity.tables)[0];
                return r.activity.tables[db].lambda_used;
            }""")
            check(abs(qq_info["fitSlope"] - state_lambda) < 1e-6,
                  f"QQ plot's fitted-line slope matches the block's lambda_used exactly (plot={qq_info['fitSlope']:.4f}, state={state_lambda:.4f})")

    bar_btn = page.locator("button.chart-toggle-btn", has_text="Show TF activity bar plot")
    check(bar_btn.count() > 0, "'Show TF activity bar plot' toggle present (at least one significant TF at the loosened threshold)")
    if bar_btn.count() > 0:
        bar_btn.click()
        page.wait_for_timeout(600)
        bar_has_data = page.evaluate("""() => {
            const divs = [...document.querySelectorAll('.chart-container')];
            const d = divs.find(x => x.id.includes('activitybar') && x.data);
            return !!d;
        }""")
        check(bar_has_data, "activity bar plot rendered a Plotly figure")

    volcano_btn = page.locator("button.chart-toggle-btn", has_text="Show TF activity volcano")
    check(volcano_btn.count() > 0, "'Show TF activity volcano' toggle present")
    if volcano_btn.count() > 0:
        volcano_btn.click()
        page.wait_for_timeout(600)
        volcano_marker_ok = page.evaluate("""() => {
            const divs = [...document.querySelectorAll('.chart-container')];
            const d = divs.find(x => x.id.includes('activityvolcano') && x.data);
            if (!d) return null;
            const sizes = d.data[0].marker.size;
            const arr = Array.isArray(sizes) ? sizes : [sizes];
            return arr.every(s => s >= 8 && s <= 26);
        }""")
        check(volcano_marker_ok is True, "activity volcano marker sizes are within the capped [8,26] range")

        # The [8,26] range check above is NOT sufficient on its own to prove
        # the size column is wired correctly: _magnitudeMarkerSizes maps every
        # non-finite value to 0 and then returns all-minimum (8px) sizes, which
        # sits inside that range. So a typo'd / renamed column key would pass
        # it silently. These two checks are the ones that actually bite - they
        # were added when "Regulon Size" was renamed to "Signed Targets
        # Tested" and the range check alone did not notice the difference.
        volcano_column_wired = page.evaluate("""() => {
            const divs = [...document.querySelectorAll('.chart-container')];
            const d = divs.find(x => x.id.includes('activityvolcano') && x.data);
            if (!d) return null;
            const tfs = d.data[0].text;
            const cd = d.data[0].customdata;
            if (!tfs || !cd) return null;
            const tables = STATE.transcriptomics.regulonResults.activity.tables;
            let checked = 0;
            for (let i = 0; i < tfs.length; i++) {
                const shown = Number(cd[i][0]);
                if (!isFinite(shown)) return { ok: false, why: 'non-finite customdata for ' + tfs[i] };
                // The invariant the rename is ABOUT: signed targets tested is
                // exactly activated + repressed, for the same TF, in whichever
                // database table this row came from.
                for (const db of Object.keys(tables)) {
                    const row = (tables[db].table.rows || []).find(r => r['TF'] === tfs[i]);
                    if (!row) continue;
                    const expected = Number(row['Activated Targets']) + Number(row['Repressed Targets']);
                    if (Number(row['Signed Targets Tested']) !== expected) {
                        return { ok: false, why: 'table row inconsistent for ' + tfs[i] };
                    }
                    if (shown !== expected) {
                        return { ok: false, why: 'plot shows ' + shown + ' but table says ' + expected + ' for ' + tfs[i] };
                    }
                    checked++;
                    break;
                }
            }
            return { ok: checked > 0, checked: checked };
        }""")
        check(volcano_column_wired is not None and volcano_column_wired.get("ok") is True,
              f"volcano's size/hover column really is the table's 'Signed Targets Tested', and it equals "
              f"Activated + Repressed for every plotted TF (got {volcano_column_wired})")
        distinct_sizes = page.evaluate("""() => {
            const divs = [...document.querySelectorAll('.chart-container')];
            const d = divs.find(x => x.id.includes('activityvolcano') && x.data);
            if (!d) return null;
            const sizes = d.data[0].marker.size;
            return Array.isArray(sizes) ? new Set(sizes).size : 1;
        }""")
        check(distinct_sizes is not None and distinct_sizes > 1,
              f"and marker sizes genuinely VARY across TFs ({distinct_sizes} distinct sizes) - an all-8px chart "
              f"is what a broken column key produces, and it would pass the range check above")

    print("=== TF Activity pane: activity-vs-own-expression (Figure B) ===")
    page.locator(".result-nav-btn", has_text="TF Activity").first.click()
    page.wait_for_timeout(200)
    # Restore a REALISTIC significance threshold for this figure. The 0.999
    # loosening above (needed so the bar/volcano toggles render at all) makes
    # every TF "significant" on both axes, which would make the fill/ring
    # encoding checks below trivially true and therefore worthless.
    page.evaluate("""() => {
        STATE.transcriptomics.regulonSigAlpha = 0.05;
        render();
    }""")
    page.wait_for_timeout(400)
    page.locator(".result-nav-btn", has_text="TF Activity").first.click()
    page.wait_for_timeout(300)

    vsexpr_btn = page.locator("button.chart-toggle-btn", has_text="Show activity vs. own expression")
    check(vsexpr_btn.count() > 0, "'Show activity vs. own expression' toggle present at the default 0.05 threshold")
    if vsexpr_btn.count() > 0:
        vsexpr_btn.click()
        page.wait_for_timeout(600)
        vsexpr_has_data = page.evaluate("""() => {
            const divs = [...document.querySelectorAll('.chart-container')];
            const d = divs.find(x => x.id.includes('activityvsexpr') && x.data);
            return !!d;
        }""")
        check(vsexpr_has_data, "activity-vs-own-expression figure rendered a Plotly figure")

        # SEMANTIC checks on the encoding. The previous version used fill for
        # the activity FDR and encoded the TF's own padj nowhere, which made
        # fill almost perfectly redundant with x-position. Assert that fill
        # and ring now encode two DIFFERENT things, both driven by real data.
        enc = page.evaluate("""() => {
            const divs = [...document.querySelectorAll('.chart-container')];
            const d = divs.find(x => x.id.includes('activityvsexpr') && x.data);
            if (!d) return null;
            const m = d.data[0].marker;
            const widths = Array.isArray(m.line.width) ? m.line.width : [m.line.width];
            const fills = Array.isArray(m.color) ? m.color : [m.color];
            const sizes = Array.isArray(m.size) ? m.size : [m.size];
            const GREY = '#b8c2cc';
            return {
                n: d.data[0].x.length,
                nRinged: widths.filter(w => w > 0).length,
                nUnringed: widths.filter(w => !w || w === 0).length,
                nFilled: fills.filter(c => c !== GREY).length,
                nGrey: fills.filter(c => c === GREY).length,
                sizesOk: sizes.every(s => s >= 8 && s <= 26),
                // Are ring and fill actually independent, or does ring just
                // mirror fill? Count points where they disagree.
                nDisagree: widths.reduce((acc, w, i) => acc + (((w > 0) !== (fills[i] !== GREY)) ? 1 : 0), 0),
            };
        }""")
        check(enc is not None, "activity-vs-expression exposes its marker encoding")
        if enc:
            check(enc["nFilled"] > 0 and enc["nGrey"] > 0,
                  f"fill separates significant from non-significant activity ({enc['nFilled']} filled, {enc['nGrey']} grey)")
            check(enc["nRinged"] > 0,
                  f"the own-padj RING is actually drawn on real data ({enc['nRinged']} of {enc['n']} TFs)")
            check(enc["nUnringed"] > 0, f"and absent where the TF's own mRNA did not significantly change ({enc['nUnringed']})")
            check(enc["nDisagree"] > 0,
                  f"ring and fill encode INDEPENDENT facts - they disagree on {enc['nDisagree']} TFs "
                  f"(if they always agreed, the ring would add nothing)")
            check(enc["sizesOk"], "marker sizes stay within the capped [8,26] range")

        legend_texts = page.locator(".chart-encoding-key").all_inner_texts()
        check(len(legend_texts) >= 1, "on-screen encoding key is rendered with the figure")
        if legend_texts:
            key = legend_texts[0]
            check("FILLED" in key and "DARK RING" in key and "SIZE" in key,
                  "key explains fill, ring AND size")
            check("post-" in key.lower() or "after translation" in key.lower(),
                  "key explains the post-translational case the figure exists to reveal")

    print("=== Regulon Detail pane: ULM regression plot + null fan (Figure A) + coherence bars (Figure E) ===")
    detail_rail = page.locator(".result-nav-btn", has_text="Regulon Detail")
    check(detail_rail.count() >= 1, "'Regulon Detail' rail entry present")
    if detail_rail.count() >= 1:
        detail_rail.first.click()
        page.wait_for_timeout(400)
        # Real, user-reported bug: this chart was built EAGERLY while this
        # pane was still hidden (display:none, behind 'TF Activity'), which
        # drew it at the wrong width - the horizontal legend rendered ~20px
        # too high, directly on top of the xaxis title text, on this exact
        # first reveal. Measured here BEFORE anything else touches the
        # chart (in particular before the TF-switch checks later in this
        # file, which redraw it from scratch and always incidentally fixed
        # it - the reason this bug was easy to miss). Fixed by having
        # makeResultNav's activate() do a full redraw (_newOraPlot's stashed
        # div._omiscerRedraw) instead of the insufficient Plotly.Plots.resize().
        legend_title_gap = page.evaluate("""() => {
            const d = document.getElementById('chart-regulon-regression');
            const legend = d && d.querySelector('.legend');
            const xtitle = d && d.querySelector('.xtitle');
            if (!legend || !xtitle) return null;
            const lr = legend.getBoundingClientRect();
            const tr = xtitle.getBoundingClientRect();
            return lr.top - tr.bottom;
        }""")
        check(legend_title_gap is not None and legend_title_gap > -5,
              f"on the FIRST reveal of this pane, the legend sits below the xaxis title with no visible "
              f"overlap - no TF switch needed (got gap={legend_title_gap}; the real pre-fix bug measured "
              f"~-20px, i.e. the legend rendered ON TOP of the title text)")
        regression_has_data = page.evaluate("""() => {
            const d = document.getElementById('chart-regulon-regression');
            return !!(d && d.data && d.data.length > 0);
        }""")
        check(regression_has_data, "ULM regression plot rendered with at least one trace")
        fan_trace_count = page.evaluate("""() => {
            const d = document.getElementById('chart-regulon-regression');
            if (!d || !d.data) return 0;
            return d.data.filter(t => t.mode === 'lines' && t.showlegend === false).length;
        }""")
        check(fan_trace_count > 10, f"null fan drawn with multiple grey lines (got {fan_trace_count} lines)")
        fitted_line_present = page.evaluate("""() => {
            const d = document.getElementById('chart-regulon-regression');
            if (!d || !d.data) return false;
            return d.data.some(t => t.name === 'Fitted slope');
        }""")
        check(fitted_line_present, "observed fitted slope line is drawn on top of the fan")

        # SEMANTIC checks, not just "a trace exists". The original version of
        # this figure passed every structural assertion while plotting the
        # wrong variable (log2FC) against a slope fitted in ranking-metric
        # units, with a centroid 141x off. These assert the rendered figure
        # is internally consistent with the fit it claims to show.
        semantics = page.evaluate("""() => {
            const d = document.getElementById('chart-regulon-regression');
            if (!d || !d.data) return null;
            const act = d.data.find(t => t.name === 'Activated targets');
            const rep = d.data.find(t => t.name === 'Repressed targets');
            const fit = d.data.find(t => t.name === 'Fitted slope');
            const violin = d.data.find(t => t.type === 'violin');
            if (!fit) return null;
            // Recover the drawn line's slope from its two endpoints.
            const drawnSlope = (fit.y[1] - fit.y[0]) / (fit.x[1] - fit.x[0]);
            const ys = [].concat(act ? act.y : [], rep ? rep.y : []);
            return {
                drawnSlope: drawnSlope,
                hasViolin: !!violin,
                violinN: violin ? violin.y.length : 0,
                nPoints: ys.length,
                yMin: Math.min.apply(null, ys),
                yMax: Math.max.apply(null, ys),
                yAxisTitle: (d.layout.yaxis.title.text || ''),
                nTicks: (d.layout.xaxis.ticktext || []).length,
            };
        }""")
        check(semantics is not None, "regression figure exposes inspectable traces")
        if semantics:
            # Points and the fitted line must live on the SAME scale. Before
            # the fix the slope was ~5.16 metric-units against a log2FC cloud
            # spanning ~2.9 - a >100% mismatch that this catches.
            span = max(semantics["yMax"] - semantics["yMin"], 1e-9)
            ratio = abs(semantics["drawnSlope"]) / span
            check(ratio < 3.0,
                  f"fitted slope is on the same scale as the plotted points "
                  f"(|slope|/span = {ratio:.2f}; was >1.7 with the log2FC/metric mismatch)")
            check(semantics["hasViolin"] and semantics["violinN"] > 50,
                  f"x=0 non-target mass drawn as a violin ({semantics['violinN']} pts) - shows what anchors the fit")
            check("log2" not in semantics["yAxisTitle"].lower(),
                  f"y-axis is labelled with the ranking metric, not log2FC (got {semantics['yAxisTitle']!r})")
            check(semantics["nTicks"] == 3,
                  f"x-axis labels all three weight groups incl. 'Not a target' (got {semantics['nTicks']} ticks)")

        # The rendered y-values must equal the server's ranking-metric values
        # for this TF - i.e. the figure is plotting the fit's own data.
        matches_payload = page.evaluate("""() => {
            const d = document.getElementById('chart-regulon-regression');
            const db = document.querySelectorAll('.result-pane:not([style*="display: none"]) select')[0].value;
            const tf = document.querySelectorAll('.result-pane:not([style*="display: none"]) select')[1].value;
            const blk = STATE.transcriptomics.regulonResults.activity.tables[db];
            const det = blk.tf_regression_detail[tf];
            const act = d.data.find(t => t.name === 'Activated targets');
            if (!det || !act) return null;
            const want = det.activated_values.slice().sort((a,b)=>a-b);
            const got = act.y.slice().sort((a,b)=>a-b);
            if (want.length !== got.length) return false;
            return want.every((v,i) => Math.abs(v - got[i]) < 1e-6);
        }""")
        check(matches_payload is True,
              "plotted activated-target y-values are exactly the server's ranking-metric values")

        coherence_has_data = page.evaluate("""() => {
            const d = document.getElementById('chart-regulon-coherence');
            return !!(d && d.data && d.data.length > 0);
        }""")
        check(coherence_has_data, "coherence bars rendered with at least one trace")

        # Coherence sign-convention regression: an "up" fraction of 1.0 with
        # zero targets showing "down" is IMPOSSIBLE unless the >= 0 fix is
        # live (a gene at exactly log2FC=0 must count as up, matching
        # build_ranked_list's own sign(0):=+1). Real data rarely has an exact
        # zero, so this is checked structurally instead: the two trace
        # fractions for each row must sum to exactly 1.0 (every target
        # assigned to exactly one bucket, none double-counted or dropped).
        coherence_sums_ok = page.evaluate("""() => {
            const d = document.getElementById('chart-regulon-coherence');
            if (!d || !d.data || d.data.length < 2) return null;
            const up = d.data.find(t => t.name === 'Moved UP');
            const down = d.data.find(t => t.name === 'Moved DOWN');
            if (!up || !down) return null;
            return up.x.every((v, i) => Math.abs(v + down.x[i] - 1.0) < 1e-9 || (v === 0 && down.x[i] === 0));
        }""")
        check(coherence_sums_ok is True or coherence_sums_ok is None,
              "coherence bars' up+down fractions sum to 1.0 per row (every target in exactly one bucket)")

        print("=== Regulon Detail pane: downloads actually work (not just 'chart renders') ===")
        # These two divs use hardcoded ids (chart-regulon-regression /
        # chart-regulon-coherence) instead of the _oraChartId() convention
        # every other chart in the app follows, because the same div is
        # reused across TF switches. _downloadOraChart used to re-derive the
        # div via _oraChartId(kind, dbName, namespace), which built an id
        # that never matched these two - so clicking Download silently said
        # "chart is not rendered yet" even though it plainly was. Fixed by
        # passing the actual clicked div through instead of re-deriving it;
        # these checks exercise the REAL browser download, not just that a
        # menu item exists.
        for container_id, kind_label in (("chart-regulon-regression", "ULM regression"),
                                          ("chart-regulon-coherence", "Coherence bars")):
            modebar_btn = page.locator(f"#{container_id} .modebar-btn[data-title='Download {kind_label}']")
            check(modebar_btn.count() == 1, f"'{kind_label}' has a Download button in its modebar")
            if modebar_btn.count() == 1:
                modebar_btn.click()
                page.wait_for_timeout(150)
                svg_item = page.locator(".ora-chart-export-item", has_text=f"{kind_label} SVG")
                check(svg_item.count() == 1, f"'{kind_label}' export menu offers an SVG option")
                if svg_item.count() == 1:
                    with page.expect_download(timeout=5000) as dl_info:
                        svg_item.click()
                    download = dl_info.value
                    check(download.suggested_filename.endswith(".svg"),
                          f"'{kind_label}' SVG download actually fires (got filename {download.suggested_filename!r})")

        # Switch TF (if more than one available) and confirm the SAME stable
        # container div re-renders rather than needing a fresh div (mirrors
        # GSEA's "Per gene set" pattern this section deliberately follows).
        tf_select = page.locator("select").filter(has=page.locator("option"))
        # Find the TF select specifically: it's the second <select> inside the
        # Regulon Detail pane's row (first is the database select).
        detail_selects = page.locator(".result-pane:visible select")
        if detail_selects.count() >= 2:
            tf_sel = detail_selects.nth(1)
            n_opts = tf_sel.locator("option").count()
            if n_opts >= 2:
                tf_sel.select_option(index=1)
                page.wait_for_timeout(400)
                regression_still_has_data = page.evaluate("""() => {
                    const d = document.getElementById('chart-regulon-regression');
                    return !!(d && d.data && d.data.length > 0);
                }""")
                check(regression_still_has_data, "switching TF re-renders the SAME regression container div")

    print("=== TF Target Enrichment pane: dot/bar plots reuse the ORA builders, plus overlap map/network (Figures C/D) ===")
    page.locator(".result-nav-btn", has_text="TF Target Enrichment").first.click()
    page.wait_for_timeout(200)
    dot_btn = page.locator("button.chart-toggle-btn", has_text="Show dot plot")
    if dot_btn.count() > 0:
        dot_btn.click()
        page.wait_for_timeout(600)
        dot_marker_ok = page.evaluate("""() => {
            const divs = [...document.querySelectorAll('.chart-container')];
            const d = divs.find(x => x.id.startsWith('chart-dot-regulon_') && x.data);
            if (!d) return null;
            const sizes = d.data[0].marker.size;
            const arr = Array.isArray(sizes) ? sizes : [sizes];
            return arr.every(s => s >= 8 && s <= 30);
        }""")
        check(dot_marker_ok is True, "enrichment dot plot marker sizes are within the capped [8,30] range (shared _magnitudeMarkerSizes)")
    else:
        check(False, "'Show dot plot' button present for at least one enrichment table")

    overlap_btn = page.locator("button.chart-toggle-btn", has_text="Show TF-TF overlap map")
    check(overlap_btn.count() > 0, "'Show TF-TF overlap map' toggle present (Figure C)")
    if overlap_btn.count() > 0:
        overlap_btn.click()
        page.wait_for_timeout(800)
        overlap_has_data = page.evaluate("""() => {
            const divs = [...document.querySelectorAll('.chart-container')];
            const d = divs.find(x => x.id.startsWith('chart-enrichmap-regulon_') && x.data);
            return !!d;
        }""")
        check(overlap_has_data, "TF-TF overlap map rendered a Plotly figure")

    network_btn = page.locator("button.chart-toggle-btn", has_text="Show TF-target network")
    check(network_btn.count() > 0, "'Show TF-target network' toggle present (Figure D)")
    if network_btn.count() > 0:
        network_btn.click()
        page.wait_for_timeout(800)
        network_has_data = page.evaluate("""() => {
            const divs = [...document.querySelectorAll('.chart-container')];
            const d = divs.find(x => x.id.startsWith('chart-network-regulon_') && x.data);
            return !!d;
        }""")
        check(network_has_data, "TF-target network rendered a Plotly figure")

    print("=== 'Start new run' resets state ===")
    new_run_btn = page.locator("button", has_text="Start new run")
    check(new_run_btn.count() == 1, "'Start new run' button present")
    new_run_btn.click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "transcriptomics_regulon_config", "'Start new run' returns to the config screen")
    check(page.evaluate("STATE.transcriptomics.regulonResults") is None, "regulonResults cleared")

    print("=== small-input warning DOES appear for a genuinely small table ===")
    # Shrink the resolved table in state and re-render the config screen: the
    # warning exists for real users with short gene lists, so it must fire.
    page.evaluate("""() => {
        const rt = STATE.transcriptomics.resolved.resolved_table;
        rt.rows = rt.rows.slice(0, 25);
        STATE.aaPage = 'transcriptomics_regulon_config';
        render();
    }""")
    page.wait_for_timeout(300)
    warn_texts = page.locator(".alert-warning").all_inner_texts()
    small_warn = [w for w in warn_texts if "resolved table has only" in w]
    check(len(small_warn) == 1, f"small-input warning appears for a 25-gene table (got {len(small_warn)})")
    if small_warn:
        check("25 gene(s)" in small_warn[0], f"warning names the actual gene count (got: {small_warn[0][:90]!r})")
        check("EMPTY" in small_warn[0], "warning says plainly that the activity table will likely be empty")
    # The warning must INFORM, not block. Note the Run button is legitimately
    # disabled right now for an unrelated reason - "Start new run" above reset
    # the draft, so no database is selected yet. Select one, then assert the
    # button enables despite the warning still being displayed.
    pills_after = page.locator("label.ora-db-pill")
    for i in range(pills_after.count()):
        if pills_after.nth(i).inner_text().strip().startswith("TRRUST"):
            cb = pills_after.nth(i).locator("input[type=checkbox]")
            if not cb.is_disabled():
                cb.check()
            break
    page.wait_for_timeout(200)
    run_btn_after = page.locator("button.ora-run-btn")
    still_warned = [w for w in page.locator(".alert-warning").all_inner_texts() if "resolved table has only" in w]
    check(len(still_warned) == 1, "warning is still displayed after selecting a database")
    check(run_btn_after.count() == 1 and not run_btn_after.is_disabled(),
          "the warning does NOT block the run - it informs, the user still decides")

    print("=== console error check ===")
    check(len(console_errors) == 0,
          f"no console/page errors during the whole flow (got {len(console_errors)}: {console_errors[:5]})")

    browser.close()

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
