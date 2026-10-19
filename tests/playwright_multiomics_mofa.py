"""Real-browser test of Multi-omics MOFA / Factor-based Integration (static/app.js).

Drives the whole flow against a live server on port 8009 with real network access (not mocked):
Multi-omics tile screen -> the MOFA box -> choose the omics -> one matrix upload each (the three 36-sample MOFA examples) -> resolve
and match -> Continue -> setup -> run -> every result pane -> a real pathway enrichment.

What is asserted, beyond "it renders":
  - the tile is live; the screen runs in MOFA mode (no correlation setup is offered, a Continue button is);
  - the setup offers the Group column of the metabolite file as the sample groups;
  - the three planted factors are recovered (|r| > 0.95 with the truth), shown by the Overview table;
  - every figure carries the figure tools (zoom, pan, in, out, autoscale, reset), a Size menu and Download figure; a bigger Size gives
    a wider figure; the figure's data are what the tables say;
  - the group table finds the treatment factor only; the factor-by-group figure has one box per group;
  - the weights figure shows 25 bars, equal to the table below it;
  - the pathway enrichment of the treatment factor puts the planted Reactome pathways (citric acid cycle / glycolysis) among the best;
  - the downloads are real CSVs with the right number of rows;
  - no page error during any of it.
Set MOFA_SHOT_DIR to save a screenshot of every pane.
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


SHOTS = os.environ.get("MOFA_SHOT_DIR")


def open_pane(page, title):
    page.locator("button.result-nav-btn", has_text=title).first.click()
    page.wait_for_timeout(900)
    if SHOTS:
        page.screenshot(path=os.path.join(SHOTS, title.replace(" ", "_") + ".png"), full_page=True)


FIG = """(id) => { const d = document.getElementById(id); return d && d.data ? { n: d.data.length, w: d._fullLayout.width, h: d._fullLayout.height,
    types: d.data.map(t => t.type), names: d.data.map(t => t.name), n0: (d.data[0].x || d.data[0].y || []).length, hasModebar: !!d.querySelector('.modebar') } : null; }"""


def tools_ok(page, fig_id):
    """The visible tools row and Download figure that sit with a figure in the active pane."""
    return page.evaluate("""(id) => {
        const d = document.getElementById(id);
        if (!d) return null;
        const pane = d.closest('.result-nav-section') || d.closest('section') || document;
        const q = (c) => [...pane.querySelectorAll('.' + c)].some(b => b.offsetParent !== null);
        return { zoom: q('fig-tool-zoom'), pan: q('fig-tool-pan'), zin: q('fig-tool-in'), zout: q('fig-tool-out'), auto: q('fig-tool-auto'), reset: q('fig-tool-reset'),
                 download: q('multiomics-mofa-download'), size: [...pane.querySelectorAll('select.multiomics-mofa-size')].some(b => b.offsetParent !== null) };
    }""", fig_id)


def all_tools(t):
    return bool(t) and all(t.values())


truth = {}
with open(os.path.join(_ROOT, "tests", "fixtures", "mofa_example_truth.csv"), newline="") as fh:
    for row in csv.DictReader(fh):
        truth[row["Sample"]] = row

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

    # ----------------------------------------------------------- tile
    print("=== the MOFA box ===")
    page.evaluate("""() => { STATE.view = 'annotation_analytics'; STATE.aaPage = 'multiomics_home'; render(); }""")
    page.wait_for_timeout(200)
    tiles = page.locator(".aa-genomics-grid button.module-btn")
    check(tiles.nth(4).inner_text().strip() == "MOFA / Factor-based Integration" and not tiles.nth(4).is_disabled(), "the fifth box is live and shows only the method name")
    check(not tiles.nth(5).is_disabled(), "SENSITIVITY: DIABLO (the sixth box) is built now too - nothing is greyed out")
    tiles.nth(4).click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "multiomics_corr" and page.evaluate("STATE.multiomics.method") == "mofa", "a fresh session opens the upload screen in MOFA mode")
    body = page.inner_text("body")
    check("MOFA / Factor-based Integration" in body and "hidden factors" in body and "Correlation-based Integration" not in body.split("MOFA / Factor-based Integration", 1)[1][:300],
          "the screen is titled for MOFA, not for correlation")

    # ----------------------------------------------------------- uploads
    print("=== uploads and sample matching ===")
    pill(page, "Transcriptomics").click()
    page.wait_for_timeout(300)
    pill(page, "Proteomics").click()
    page.wait_for_timeout(500)
    pill(page, "Metabolomics").click()
    page.wait_for_timeout(500)
    check(page.locator("input[type=file]").count() == 3, "one upload box per chosen omics")
    page.evaluate("""() => { STATE.multiomics.corr.taxonId = '9606'; STATE.multiomics.corr.taxonName = 'Homo sapiens'; render(); }""")
    page.wait_for_timeout(200)
    for i, f in enumerate(["example_mofa_transcriptomics.csv", "example_mofa_proteomics.csv", "example_mofa_metabolomics.csv"]):
        page.locator("input[type=file]").nth(i).set_input_files(example(f))
        page.wait_for_timeout(1500)
    with page.expect_download() as dl:
        page.locator("button.ora-inline-example-link").first.click()
    check(dl.value.suggested_filename == "example_mofa_transcriptomics.csv", f"the example link serves the 36-sample MOFA example ({dl.value.suggested_filename})")
    check(page.evaluate("STATE.multiomics.corr.uploads.metabolomics.info.set_aside_columns") == ["Group"], "the metabolite table's Group column is set aside, not data")
    page.locator("button.multiomics-corr-resolve").click()
    wait_for(page, "!!STATE.multiomics.corr.alignment", 240, "identifiers resolved and samples matched")
    check(page.evaluate("STATE.multiomics.corr.alignment.n_shared") == 36, "36 samples are in every file")
    check(page.locator("button.multiomics-corr-run").count() == 0, "no correlation setup is offered in MOFA mode")
    cont = page.locator("button.multiomics-mofa-continue-btn")
    check(cont.count() == 1 and not cont.is_disabled(), "'Continue to MOFA setup' is offered")

    # ----------------------------------------------------------- setup
    print("=== setup ===")
    cont.click()
    wait_for(page, "!!STATE.multiomics.mofa.options", 20, "options loaded")
    page.wait_for_timeout(400)
    check(page.evaluate("STATE.aaPage") == "multiomics_mofa", "Continue opens the MOFA setup")
    check(page.locator("select.multiomics-mofa-group").evaluate("s => s.value") == "metabolomics:Group", "the Group column of the metabolite file is the default sample group")
    gtxt = page.locator("select.multiomics-mofa-group").inner_text()
    check("Control" in gtxt and "Treated" in gtxt and "36 of 36" in gtxt, f"the group source says which groups and how many samples ({gtxt[:90]!r})")
    st = page.inner_text("body")
    check("Offer MORE factors" in st and "NOT used" not in st and "ALL samples together" in st, "the setup explains that groups are not used to fit the model and to offer more factors")
    check(page.evaluate("STATE.multiomics.mofa.draft.factors") == 10 and page.evaluate("STATE.multiomics.mofa.draft.n_refits") == 3, "defaults: 10 factors, 3 refits")
    if SHOTS:
        page.screenshot(path=os.path.join(SHOTS, "setup.png"), full_page=True)
    page.locator("button.multiomics-mofa-run").click()
    wait_for(page, "!!STATE.multiomics.mofa.run.results", 240, "MOFA finished")
    page.wait_for_timeout(800)
    R = page.evaluate("STATE.multiomics.mofa.run.results")
    check(R is not None and page.evaluate("STATE.aaPage") == "multiomics_mofa_results", "the run reached the results screen")

    # ----------------------------------------------------------- overview
    print("=== overview ===")
    fv = R["factor_values"]
    Z = np.array(fv["values"])
    rec = {}
    for name in ("treatment", "cholesterol", "lipids"):
        t = np.array([float(truth[s][name]) for s in fv["samples"]])
        rec[name] = max(abs(np.corrcoef(Z[:, k], t)[0, 1]) for k in range(Z.shape[1]))
    check(all(v > 0.95 for v in rec.values()), f"SENSITIVITY: the three planted factors are recovered ({ {k: round(v, 3) for k, v in rec.items()} })")
    check(R["converged"] and R["n_factors"] >= 3, f"converged, {R['n_factors']} factors kept")
    cards = page.inner_text(".metric-cards, .metric-card-row") if page.locator(".metric-cards, .metric-card-row").count() else page.inner_text("body")
    check("36" in cards and "Transcriptomics + Proteomics + Metabolomics" in cards, "the cards name the samples and the omics")
    ov = page.inner_text("body")
    check("Stability (median |r| in the refits)" in ov and "stable" in ov and "UNSTABLE" in ov, "the Overview table shows each factor's stability, the unstable one flagged")
    check("not a pathway and not a cause" in ov, "the Overview says what the results do not tell you")
    if SHOTS:
        page.screenshot(path=os.path.join(SHOTS, "Overview.png"), full_page=True)

    # ----------------------------------------------------------- variance explained
    print("=== variance explained ===")
    open_pane(page, "Variance explained")
    f = page.evaluate(FIG, "multiomics-mofa-r2")
    check(f and f["types"] == ["heatmap"] and f["n0"] == R["n_factors"], "the variance heatmap is drawn, one column per factor")
    z = page.evaluate("document.getElementById('multiomics-mofa-r2').data[0].z")
    check(np.allclose(np.array(z), np.array(R["r2"]["values"])), "its cells are the table's percentages")
    check(all_tools(tools_ok(page, "multiomics-mofa-r2")), f"it has the figure tools, Size and Download figure ({tools_ok(page, 'multiomics-mofa-r2')})")
    f2 = page.evaluate(FIG, "multiomics-mofa-totals")
    check(f2 and f2["types"] == ["bar"] and f2["n0"] == 3 and all_tools(tools_ok(page, "multiomics-mofa-totals")), "the totals figure has three bars and the same tools")
    w100 = f["w"]
    page.locator("select.multiomics-mofa-size[data-fig='multiomics-mofa-r2']").select_option("1.5")
    page.wait_for_timeout(700)
    w150 = page.evaluate(FIG, "multiomics-mofa-r2")["w"]
    check(abs(w150 / w100 - 1.5) < 0.05, f"SENSITIVITY: Size 150% makes the figure 1.5 x wider ({w100} -> {w150})")
    page.locator("select.multiomics-mofa-size[data-fig='multiomics-mofa-r2']").select_option("1")
    page.wait_for_timeout(500)
    tot_shared = np.array(R["r2"]["values"])
    shared = [k for k in range(R["n_factors"]) if (tot_shared[:, k] > 5).sum() >= 3]
    check(len(shared) >= 1, "at least one factor is shared by all three omics (the treatment factor)")
    with page.expect_download() as dl:
        page.locator("button", has_text="Download variance explained").click()
    rows = list(csv.reader(io.StringIO(open(dl.value.path(), encoding="utf-8").read())))
    check(len(rows) == 4, "the variance CSV has a header and one row per omics")

    # ----------------------------------------------------------- factors
    print("=== factors ===")
    open_pane(page, "Factors")
    f = page.evaluate(FIG, "multiomics-mofa-scatter")
    check(f and f["n"] == 2 and sorted(f["names"]) == ["Control", "Treated"] and f["n0"] == 18, "the factor scatter has one trace per group, 18 samples each")
    check(all_tools(tools_ok(page, "multiomics-mofa-scatter")), "the scatter has the figure tools, Size and Download figure")
    page.locator("select.multiomics-mofa-fy").select_option("2")
    page.wait_for_timeout(600)
    y0 = page.evaluate("document.getElementById('multiomics-mofa-scatter').data[0].y")
    idx = [i for i, s in enumerate(fv["samples"]) if s.startswith("Control")]
    check(np.allclose(y0, [Z[i, 2] for i in idx]), "changing the vertical factor redraws it with that factor's values")
    f = page.evaluate(FIG, "multiomics-mofa-box")
    check(f and f["types"] == ["box", "box"] and all_tools(tools_ok(page, "multiomics-mofa-box")), "the factor-by-group figure has one box per group and the tools")
    gt = R["group_tests"]["rows"]
    sig = [r["Factor"] for r in gt if r["FDR"] < 0.05]
    k_treat = int(np.argmax([abs(np.corrcoef(Z[:, k], np.array([float(truth[s]["treatment"]) for s in fv["samples"]]))[0, 1]) for k in range(Z.shape[1])]))
    check(sig == [fv["names"][k_treat]], f"SENSITIVITY: only the treatment factor differs between the groups ({sig})")
    txt = page.inner_text("body")
    check("Mann-Whitney U" in txt and "rank-biserial" in txt, "the group table names its test and its effect measure")

    # ----------------------------------------------------------- weights
    print("=== weights ===")
    open_pane(page, "Weights")
    f = page.evaluate(FIG, "multiomics-mofa-weights")
    check(f and f["types"] == ["bar"] and f["n0"] == 25 and all_tools(tools_ok(page, "multiomics-mofa-weights")), "the weights figure shows 25 bars and has the tools")
    page.locator("select.multiomics-mofa-womics").select_option("proteomics")
    page.wait_for_timeout(700)
    ys = page.evaluate("document.getElementById('multiomics-mofa-weights').data[0].y")
    top = [x["label"] for x in R["top_weights"]["proteomics"][fv["names"][0]]]
    check(len(set(ys)) == len(ys) == 25, "SENSITIVITY: every bar has its own row (two features of one gene are not stacked on one row)")
    xs = page.evaluate("document.getElementById('multiomics-mofa-weights').data[0].x")
    check(all(abs(xs[i]) >= abs(xs[i + 1]) - 1e-12 for i in range(24)), "the bars run from the largest |weight| down")
    page.locator("select.multiomics-mofa-womics").select_option("transcriptomics")
    page.wait_for_timeout(700)
    ys_t = page.evaluate("document.getElementById('multiomics-mofa-weights').data[0].y")
    check(len(set(ys_t)) == 25, "...also for the transcripts, where one gene symbol is carried by two features of this example")
    page.locator("select.multiomics-mofa-womics").select_option("proteomics")
    page.wait_for_timeout(700)
    check(ys == top, "SENSITIVITY: the proteomics bars are the table's 25 largest |weights| of Factor 1, in order")
    with page.expect_download() as dl:
        page.locator("button", has_text="Download all weights").click()
    rows = list(csv.reader(io.StringIO(open(dl.value.path(), encoding="utf-8").read())))
    check(len(rows) - 1 == R["n_features"]["proteomics"], f"the weights CSV has every feature that entered the model ({len(rows) - 1})")

    # ----------------------------------------------------------- enrichment (real Reactome)
    print("=== pathway enrichment ===")
    open_pane(page, "Pathway enrichment")
    page.locator("select.multiomics-mofa-eomics").select_option("transcriptomics")
    page.locator("select.multiomics-mofa-efactor").select_option(str(k_treat))
    page.locator("button.multiomics-mofa-enrich-run").click()
    wait_for(page, "Object.keys(STATE.multiomics.mofa.enrich).length > 0", 240, "enrichment finished")
    page.wait_for_timeout(1200)
    E = page.evaluate("Object.values(STATE.multiomics.mofa.enrich)[0]")
    names = [r["Term Name"] for r in E["table"]["rows"][:8]]
    print("   top pathways:", names)
    check(E["n_tested"] > 20 and any(("citric" in n.lower() or "glycoly" in n.lower() or "tca" in n.lower()) for n in names),
          f"SENSITIVITY: the planted pathways (citric acid cycle / glycolysis) are among the best of {E['n_tested']} tested")
    check(page.locator("button", has_text="Download this enrichment").count() >= 1, "the enrichment can be downloaded")
    f = page.evaluate(FIG, "multiomics-mofa-enrich-dots")
    check(f and f["types"] == ["scatter"] and all_tools(tools_ok(page, "multiomics-mofa-enrich-dots")), "the enrichment dot plot is drawn with the tools, Size and Download figure")
    if SHOTS:
        page.screenshot(path=os.path.join(SHOTS, "enrichment.png"), full_page=True)
    page.locator("select.multiomics-mofa-emode").select_option("absolute")
    page.wait_for_timeout(400)
    check(page.locator("button", has_text="Download this enrichment").count() == 0, "a different test has no result yet (nothing from the signed run is shown for it)")
    page.locator("button.multiomics-mofa-enrich-run").click()
    wait_for(page, "Object.keys(STATE.multiomics.mofa.enrich).length > 1", 240, "absolute enrichment finished")
    page.wait_for_timeout(1200)
    E2 = page.evaluate("Object.values(STATE.multiomics.mofa.enrich).find(e => e.mode === 'absolute')")
    check(E2 is not None and E2["mode"] == "absolute" and E2["n_tested"] > 20, "SENSITIVITY: the chosen test (absolute weights) is the one that was run")
    title = page.evaluate("document.getElementById('multiomics-mofa-enrich-dots')._fullLayout.xaxis.title.text")
    check("|weight|" in title, f"its dot plot says what the AUC is about ({title[:60]})")
    page.locator("select.multiomics-mofa-emode").select_option("signed")
    page.wait_for_timeout(500)
    check(page.locator("button", has_text="Download this enrichment").count() >= 1, "switching back to the signed test shows its earlier result")

    # ----------------------------------------------------------- data used, parameters, downloads
    print("=== data used, parameters, downloads ===")
    open_pane(page, "Data used")
    check("Too many missing values" in page.inner_text("body"), "the Data used pane reports the features dropped")
    open_pane(page, "Parameters used")
    check("mofapy2" in page.inner_text("body") and "Random seed: 42" in page.inner_text("body"), "the parameters name the model and the seed")
    open_pane(page, "Factors")
    with page.expect_download() as dl:
        page.locator("button", has_text="Download factor values per sample").click()
    rows = list(csv.reader(io.StringIO(open(dl.value.path(), encoding="utf-8").read())))
    check(len(rows) == 37, f"the factor-values CSV has one row per sample ({len(rows) - 1})")

    # ----------------------------------------------------------- coming back
    print("=== coming back and changing the tables ===")
    page.evaluate("""() => { STATE.aaPage = 'multiomics_home'; render(); }""")
    page.wait_for_timeout(200)
    page.locator(".aa-genomics-grid button.module-btn").nth(4).click()
    page.wait_for_timeout(400)
    check(page.evaluate("STATE.aaPage") == "multiomics_corr" and page.evaluate("!!STATE.multiomics.corr.alignment") and page.locator("button.multiomics-mofa-continue-btn").count() == 1,
          "choosing the MOFA box again keeps the tables that are already matched")
    check(page.evaluate("!!STATE.multiomics.mofa.run.results"), "...and the MOFA result is still there")
    page.locator("input[type=file]").nth(0).set_input_files(example("example_mofa_transcriptomics.csv"))
    page.wait_for_timeout(1500)
    check(page.evaluate("STATE.multiomics.corr.alignment") is None and page.evaluate("STATE.multiomics.mofa.run.results") is None
          and page.evaluate("STATE.multiomics.mofa.options") is None, "SENSITIVITY: uploading a table again clears the match AND the MOFA result")

    check(len(errors) == 0, f"no page error during the whole flow (got {errors[:3]})")
    browser.close()

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
