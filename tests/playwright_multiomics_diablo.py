"""Real-browser test of Multi-omics DIABLO / Multiblock Integration (static/app.js).

Drives the whole flow against a live server on port 8009 with real network access: Multi-omics tile screen -> the DIABLO box
-> choose the omics -> one matrix upload each (the 36-sample MOFA examples) -> resolve and match -> Continue -> setup (group
REQUIRED) -> run (component search, keepX search, cross-validated error, permutation test) -> every result pane -> a real
pathway enrichment.
"""
import csv
import io
import os

import numpy as np
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


def wait_for(page, expr, seconds, label):
    for _ in range(seconds):
        if page.evaluate(expr):
            return True
        page.wait_for_timeout(1000)
    check(False, f"timed out waiting for: {label}")
    return False


def example(name):
    return os.path.join(_ROOT, "Data3", name)


def pill(page, label):
    return page.locator("label.ora-db-pill", has_text=label).first.locator("input")


SHOTS = os.environ.get("DIABLO_SHOT_DIR")


def open_pane(page, title):
    page.locator("button.result-nav-btn", has_text=title).first.click()
    page.wait_for_timeout(900)
    if SHOTS:
        page.screenshot(path=os.path.join(SHOTS, title.replace(" ", "_") + ".png"), full_page=True)


FIG = """(id) => { const d = document.getElementById(id); return d && d.data ? { n: d.data.length, w: d._fullLayout.width, h: d._fullLayout.height,
    types: d.data.map(t => t.type) } : null; }"""


def tools_ok(page, fig_id):
    return page.evaluate("""(id) => {
        const d = document.getElementById(id);
        if (!d) return null;
        const pane = d.closest('.result-nav-section') || d.closest('section') || document;
        const q = (c) => [...pane.querySelectorAll('.' + c)].some(b => b.offsetParent !== null);
        return { zoom: q('fig-tool-zoom'), pan: q('fig-tool-pan'), zin: q('fig-tool-in'), zout: q('fig-tool-out'), auto: q('fig-tool-auto'), reset: q('fig-tool-reset'),
                 download: q('multiomics-diablo-download'), size: [...pane.querySelectorAll('select.multiomics-diablo-size')].some(b => b.offsetParent !== null) };
    }""", fig_id)


def all_tools(t):
    return bool(t) and all(t.values())


with sync_playwright() as p:
    browser = p.chromium.launch()
    context = browser.new_context(accept_downloads=True, viewport={"width": 1400, "height": 1000})
    page = context.new_page()
    errors = []
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    page.on("console", lambda msg: errors.append("console: " + msg.text) if msg.type == "error" else None)
    page.goto(BASE, wait_until="networkidle")
    page.wait_for_timeout(300)
    check(len(errors) == 0, f"app.js parsed with no top-level error (got {errors})")

    print("=== the DIABLO box ===")
    page.evaluate("""() => { STATE.view = 'annotation_analytics'; STATE.aaPage = 'multiomics_home'; render(); }""")
    page.wait_for_timeout(200)
    tiles = page.locator(".aa-genomics-grid button.module-btn")
    check(tiles.count() == 6 and tiles.nth(5).inner_text().strip() == "DIABLO / Multiblock Integration" and not tiles.nth(5).is_disabled(),
          "the sixth (last) box is live and shows only the method name")
    tiles.nth(5).click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "multiomics_corr" and page.evaluate("STATE.multiomics.method") == "diablo", "a fresh session opens the upload screen in DIABLO mode")
    body = page.inner_text("body")
    check("DIABLO / Multiblock Integration" in body and "SUPERVISED" in body, "the screen is titled for DIABLO and says it is supervised")

    print("=== uploads and sample matching ===")
    pill(page, "Transcriptomics").click()
    page.wait_for_timeout(300)
    pill(page, "Proteomics").click()
    page.wait_for_timeout(500)
    pill(page, "Metabolomics").click()
    page.wait_for_timeout(500)
    page.evaluate("""() => { STATE.multiomics.corr.taxonId = '9606'; STATE.multiomics.corr.taxonName = 'Homo sapiens'; render(); }""")
    page.wait_for_timeout(200)
    for i, f in enumerate(["example_mofa_transcriptomics.csv", "example_mofa_proteomics.csv", "example_mofa_metabolomics.csv"]):
        page.locator("input[type=file]").nth(i).set_input_files(example(f))
        page.wait_for_timeout(1500)
    page.locator("button.multiomics-corr-resolve").click()
    wait_for(page, "!!STATE.multiomics.corr.alignment", 240, "identifiers resolved and samples matched")
    check(page.evaluate("STATE.multiomics.corr.alignment.n_shared") == 36, "36 samples are in every file")
    cont = page.locator("button.multiomics-diablo-continue-btn")
    check(cont.count() == 1 and not cont.is_disabled(), "'Continue to DIABLO setup' is offered")

    print("=== setup ===")
    cont.click()
    wait_for(page, "!!STATE.multiomics.diablo.options", 20, "options loaded")
    page.wait_for_timeout(400)
    check(page.evaluate("STATE.aaPage") == "multiomics_diablo", "Continue opens the DIABLO setup")
    check(page.locator("select.multiomics-diablo-group").evaluate("s => s.value") == "metabolomics:Group", "the Group column of the metabolite file is the default sample group")
    check(not page.locator("button.multiomics-diablo-run").is_disabled(), "Run is enabled once a group is chosen")
    # the setup screen never lets the draft's group_column go stale: if it points at a column that is no longer offered,
    # the SELECT falls back to the first usable option and the draft is corrected to match - verified directly, since
    # trying to force an empty group_column open-loop (as a first draft of this test did) gets silently auto-corrected
    # on the very next render, which is the intended behaviour, not a gap: the button can only be disabled when NO
    # usable group source exists at all (checked server-side regardless - see test_multiomics_diablo_session.py).
    check(page.evaluate("STATE.multiomics.diablo.draft.group_column") == "metabolomics:Group", "the draft's group column always holds a valid, usable option")
    # keep the run small/fast for the browser test
    page.fill("input.multiomics-diablo-grid", "5, 10")
    page.locator("input.multiomics-diablo-grid").dispatch_event("change")
    page.wait_for_timeout(200)
    check(page.evaluate("STATE.multiomics.diablo.draft.keepX_grid") == [5, 10], "the feature-count grid field updates the draft")
    if SHOTS:
        page.screenshot(path=os.path.join(SHOTS, "setup.png"), full_page=True)
    page.locator("button.multiomics-diablo-run").click()
    wait_for(page, "!!STATE.multiomics.diablo.run.results", 240, "DIABLO finished")
    page.wait_for_timeout(800)
    R = page.evaluate("STATE.multiomics.diablo.run.results")
    check(R is not None and page.evaluate("STATE.aaPage") == "multiomics_diablo_results", "the run reached the results screen")

    print("=== overview ===")
    check(R["cv"]["BER"] < 0.1, f"SENSITIVITY: the cross-validated error on the real planted treatment signal is small ({R['cv']['BER']:.3f})")
    check(R["permutation"]["p_value"] < 0.2, f"the permutation test flags it as unlikely by chance (p={R['permutation']['p_value']:.3f})")
    cards = page.inner_text("body")
    check("36" in cards and "Transcriptomics + Proteomics + Metabolomics" in cards, "the cards name the samples and the omics")
    ov = page.inner_text("body")
    check("Permutation p-value" in ov and "shuffled" in ov, "the Overview explains the permutation test")
    check("not necessarily the biologically important one" in ov, "the Overview says what the results do not tell you")
    if SHOTS:
        page.screenshot(path=os.path.join(SHOTS, "Overview.png"), full_page=True)

    print("=== components ===")
    open_pane(page, "Components")
    n_groups = len(set(R["groups"].values()))
    if R["n_components"] >= 2:
        f = page.evaluate(FIG, "multiomics-diablo-scatter")
        check(f and "scatter" in f["types"] and f["n"] == len(R["omics"]) * n_groups and all_tools(tools_ok(page, "multiomics-diablo-scatter")),
              f"SENSITIVITY: the component scatter has one trace per omics x group ({f['n'] if f else None} of {len(R['omics']) * n_groups} expected)")
        total_pts = page.evaluate("document.getElementById('multiomics-diablo-scatter').data.reduce((s, t) => s + (t.x || []).length, 0)")
        check(total_pts == len(R["omics"]) * len(R["samples"]), f"SENSITIVITY: each omics' points are split ACROSS its group traces, not duplicated into every one ({total_pts} of {len(R['omics']) * len(R['samples'])} expected)")
    f2 = page.evaluate(FIG, "multiomics-diablo-box")
    check(f2 and all(t == "box" for t in f2["types"]) and f2["n"] == len(R["omics"]) * n_groups and all_tools(tools_ok(page, "multiomics-diablo-box")),
          f"SENSITIVITY: the per-group box figure has one box per omics x group ({f2['n'] if f2 else None} of {len(R['omics']) * n_groups} expected), drawn with the tools")
    with page.expect_download() as dl:
        page.locator("button", has_text="Download component values").click()
    rows = list(csv.reader(io.StringIO(open(dl.value.path(), encoding="utf-8").read())))
    check(len(rows) == 37, f"the component-values CSV has one row per sample ({len(rows) - 1})")

    print("=== selected features ===")
    open_pane(page, "Selected features")
    f = page.evaluate(FIG, "multiomics-diablo-weights")
    check(f and f["types"] == ["bar"] and all_tools(tools_ok(page, "multiomics-diablo-weights")), "the loadings figure is drawn with the tools")
    # a duplicate LABEL (two features sharing one gene symbol - isoforms, protein groups) must get its feature id appended
    # so both stay distinguishable, as MOFA's weights chart does; recompute the expected row labels independently, for
    # whichever omics shows a genuine duplicate in this run (checked across all of them, not assumed to be one specific omics)
    mismatch = None
    for o in R["omics"]:
        page.locator("select.multiomics-diablo-womics").select_option(o)
        page.wait_for_timeout(700)
        ys = page.evaluate("document.getElementById('multiomics-diablo-weights').data[0].y")
        items = [x for x in R["top_weights"][o][R["component_names"][0]] if x["weight"] != 0][:40]
        seen = {}
        for x in items:
            seen[x["label"]] = seen.get(x["label"], 0) + 1
        expected = [(x["label"] + " [" + x["feature"] + "]") if seen[x["label"]] > 1 and x["feature"] != x["label"] else x["label"] for x in items]
        if ys != expected:
            mismatch = (o, ys, expected)
    check(mismatch is None, f"SENSITIVITY: every omics' selected-feature labels match the expected dedup rule ({mismatch})")
    has_dup = any(len(set(x["label"] for x in R["top_weights"][o][R["component_names"][0]] if x["weight"] != 0)) <
                 len([x for x in R["top_weights"][o][R["component_names"][0]] if x["weight"] != 0]) for o in R["omics"])
    print("   (a genuine duplicate label was present in this run:", has_dup, ")")
    page.locator("select.multiomics-diablo-womics").select_option("metabolomics")
    page.wait_for_timeout(700)

    # whether or not THIS run's real data happens to contain a duplicate gene-symbol label, exercise the rule directly:
    # a synthetic pair of features sharing one label must both stay distinguishable (the feature id appended), not merge
    dedup_labels = page.evaluate("""() => {
        const div = document.createElement('div');
        document.body.appendChild(div);
        const fakeR = { top_weights: { fakeomics: { C1: [
            { feature: 'ENSG1', label: 'SAMEGENE', weight: 0.9 },
            { feature: 'ENSG2', label: 'SAMEGENE', weight: -0.5 },
            { feature: 'ENSG3', label: 'UNIQUEGENE', weight: 0.3 },
        ] } } };
        buildDiabloWeights(div, fakeR, 'fakeomics', 'C1', 'test-dedup', 1);
        const ys = div.data[0].y;
        div.remove();
        return ys;
    }""")
    check(dedup_labels == ["SAMEGENE [ENSG1]", "SAMEGENE [ENSG2]", "UNIQUEGENE"],
          f"SENSITIVITY: two features sharing one gene symbol BOTH get their feature id appended so they stay distinguishable, a third with a unique label does not ({dedup_labels})")

    print("=== performance ===")
    open_pane(page, "Performance")
    f = page.evaluate(FIG, "multiomics-diablo-perf")
    check(f and f["types"] == ["bar"] and all_tools(tools_ok(page, "multiomics-diablo-perf")), "the component-count search figure is drawn with the tools")
    f2 = page.evaluate(FIG, "multiomics-diablo-perm")
    check(f2 and f2["types"] == ["histogram"] and all_tools(tools_ok(page, "multiomics-diablo-perm")), "the permutation histogram is drawn with the tools")
    null_mean = page.evaluate("document.getElementById('multiomics-diablo-perm').layout.shapes[0].x0")
    check(abs(null_mean - R["cv"]["BER"]) < 1e-9, "the observed-BER marker line sits at the real cross-validated error")
    with page.expect_download() as dl:
        page.locator("button", has_text="Download the cross-validated error").click()
    rows = list(csv.reader(io.StringIO(open(dl.value.path(), encoding="utf-8").read())))
    check(len(rows) == 3, "the CV-error CSV has the BER and ER rows")

    print("=== pathway enrichment (real Reactome) ===")
    open_pane(page, "Pathway enrichment")
    page.locator("select.multiomics-diablo-eomics").select_option("transcriptomics")
    page.locator("button.multiomics-diablo-enrich-run").click()
    wait_for(page, "Object.keys(STATE.multiomics.diablo.enrich).length > 0", 240, "enrichment finished")
    page.wait_for_timeout(1200)
    E = page.evaluate("Object.values(STATE.multiomics.diablo.enrich)[0]")
    check(E["n_tested"] > 10, f"{E['n_tested']} pathways tested")
    names = [r["Term Name"] for r in E["table"]["rows"][:8]]
    print("   top pathways:", names)
    check(page.locator("button", has_text="Download this enrichment").count() >= 1, "the enrichment can be downloaded")
    f = page.evaluate(FIG, "multiomics-diablo-enrich-dots")
    check(f and f["types"] == ["scatter"] and all_tools(tools_ok(page, "multiomics-diablo-enrich-dots")), "the enrichment dot plot is drawn with the tools")
    if SHOTS:
        page.screenshot(path=os.path.join(SHOTS, "enrichment.png"), full_page=True)

    print("=== data used, parameters ===")
    open_pane(page, "Data used")
    check("Used (most variable)" in page.inner_text("body"), "the Data used pane reports the features used")
    open_pane(page, "Parameters used")
    ptxt = page.inner_text("body")
    check("mixOmics" in ptxt and "Random seed: 42" in ptxt and "Sample groups: metabolomics:Group" in ptxt, "the parameters name the model, the seed and the group source")

    print("=== coming back and changing the tables ===")
    page.evaluate("""() => { STATE.aaPage = 'multiomics_home'; render(); }""")
    page.wait_for_timeout(200)
    page.locator(".aa-genomics-grid button.module-btn").nth(5).click()
    page.wait_for_timeout(400)
    check(page.evaluate("STATE.aaPage") == "multiomics_corr" and page.evaluate("!!STATE.multiomics.corr.alignment") and page.locator("button.multiomics-diablo-continue-btn").count() == 1,
          "choosing the DIABLO box again keeps the tables that are already matched")
    check(page.evaluate("!!STATE.multiomics.diablo.run.results"), "...and the DIABLO result is still there")
    page.locator("input[type=file]").nth(0).set_input_files(example("example_mofa_transcriptomics.csv"))
    page.wait_for_timeout(1500)
    check(page.evaluate("STATE.multiomics.corr.alignment") is None and page.evaluate("STATE.multiomics.diablo.run.results") is None,
          "SENSITIVITY: uploading a table again clears the match AND the DIABLO result")

    check(len(errors) == 0, f"no page error during the whole flow (got {errors[:3]})")
    browser.close()

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
