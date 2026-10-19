"""Real-browser test of Multi-omics Network-based Integration (static/app.js).

Drives the whole flow against a live server on port 8009 with real network access (not mocked): Multi-omics tile screen ->
the Network-based Integration box -> all four omics -> one upload each -> resolve -> Continue -> setup -> run (500
permutations on the real human STRING network) -> every result pane.

Asserted, beyond "it renders":
  - the fourth box is live, the fifth is still greyed out; the shared upload screen leads here (not to Joint's setup);
  - the setup shows only what applies (the metabolite limits appear because metabolomics is chosen), and the run time hint;
  - the overview lists every omics' contribution; the hit table opens on the RECRUITED view when there are recruited proteins;
  - the single-omics tiles' own figures (network, null distribution) draw from this payload;
  - the module composition chart has one bar series per omics (when a module exists);
  - the pathway network is drawn at a readable fixed size in a scroll frame with the figure tools row: every node drawn, no two
    nodes closer than a marker, NO TWO LABELS OVERLAP (measured on the page), one chip series per omics;
  - downloads are real CSVs; no page error.
"""
import csv
import io
import os
import re

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


def open_pane(page, title):
    page.locator("button.result-nav-btn", has_text=title).first.click()
    page.wait_for_timeout(700)
    if os.environ.get("NET_SHOT_DIR"):
        page.screenshot(path=os.path.join(os.environ["NET_SHOT_DIR"], re.sub(r"\W+", "_", title) + ".png"), full_page=True)


HAS_DATA = """(id) => { const d = document.getElementById(id); return !!(d && d.data && d.data.length); }"""
MODEBAR = """(id) => { const d = document.getElementById(id); return !!(d && d.querySelector('.modebar')); }"""

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
    print("=== the Network-based Integration box ===")
    page.evaluate("""() => { STATE.view = 'annotation_analytics'; STATE.aaPage = 'multiomics_home'; render(); }""")
    page.wait_for_timeout(200)
    tiles = page.locator(".aa-genomics-grid button.module-btn")
    check(tiles.nth(3).inner_text().strip() == "Network-based Integration" and not tiles.nth(3).is_disabled(), "the fourth box is live and shows only its name")
    check(tiles.nth(5).is_disabled() and "coming soon" in tiles.nth(5).inner_text(), "the only unbuilt box (DIABLO) is still greyed out")
    tiles.nth(3).click()
    page.wait_for_timeout(300)
    check(page.evaluate("STATE.aaPage") == "multiomics" and page.evaluate("STATE.multiomics.method") == "network", "a fresh session opens the shared upload screen in network mode")
    body = page.inner_text("body")
    check("Multi-omics — Network-based Integration" in body and "centre of ALL your omics" in body, "the screen is titled for this method (not Joint Pathway Enrichment)")

    # ----------------------------------------------------------- uploads + resolve
    print("=== four omics: upload and resolve ===")
    for label in ["Genomics", "Transcriptomics", "Proteomics", "Metabolomics"]:
        if not pill(page, label).is_checked():
            pill(page, label).click()
            page.wait_for_timeout(400)
    page.evaluate("""() => { STATE.multiomics.taxonId = '9606'; STATE.multiomics.taxonName = 'Homo sapiens'; render(); }""")
    page.wait_for_timeout(200)
    for i, f in enumerate(["example_multiomics_genes.csv", "example_transcriptomics_de_full.csv", "example_proteomics_de_full.csv", "example_metabolomics_list.csv"]):
        page.locator("input[type=file]").nth(i).set_input_files(example(f))
        page.wait_for_timeout(1800)
    page.locator("button", has_text="Resolve identifiers").first.click()
    check(wait_for(page, "!!STATE.multiomics.resolved", 900, "resolve job"), "resolve completed")
    page.wait_for_timeout(300)
    cont = page.locator("button", has_text="Continue to network setup")
    check(cont.count() == 1 and page.locator("button", has_text="Continue to setup").count() == 0, "the Continue button leads to the network setup (not to Joint's setup)")
    cont.first.click()
    page.wait_for_timeout(500)
    check(page.evaluate("STATE.aaPage") == "multiomics_net", "navigated to the network setup")

    # ----------------------------------------------------------- setup
    print("=== setup ===")
    body = page.inner_text("body")
    check("Equal weight per omics" in body and "Metabolomics: a hit metabolite seeds the proteins" in body, "the setup states how each omics takes part and the equal-weight rule")
    check(page.locator("input[type=number]").count() == 12, f"the cutoffs, FDR, module, metabolite and pathway fields are there ({page.locator('input[type=number]').count()} number fields: 4 cutoffs, FDR, 2 module, 2 metabolite, 3 pathway)")
    check("Merge a pathway into a better-ranked one" in body, "the pathway-merge threshold is a setup option")
    check("Leave out a compound found in more than" in body and "Skip a reaction that lists more than" in body, "the metabolite limits appear because metabolomics is chosen")
    check(page.locator("select.multiomics-net-confidence").input_value() == "700" and page.locator("select.multiomics-net-perm").input_value() == "1000", "defaults: STRING confidence 0.7, 1000 permutations")
    page.locator("select.multiomics-net-perm").select_option("500")
    check("about 3 minute(s)" in page.inner_text("body"), "the run-time hint follows the permutation count (500 permutations: about 3 minutes)")

    # ----------------------------------------------------------- run
    print("=== run ===")
    page.locator("button.multiomics-net-run").click()
    page.wait_for_timeout(1500)
    check(page.evaluate("STATE.aaPage") == "multiomics_net_results" and "Running analysis" in page.inner_text("body"), "the results screen opens with a progress bar")
    check(wait_for(page, "!!STATE.multiomics.net.run.results", 1500, "network run"), "the run completed")
    page.wait_for_timeout(1200)
    R = page.evaluate("STATE.multiomics.net.run.results")
    rail = [t.strip() for t in page.locator("button.result-nav-btn").all_inner_texts()]
    expected = ["Overview", "Hit table", "Figures"] + (["Modules"] if R["modules"] and R["modules"]["table"]["rows"] else []) + (
        ["Module enrichment"] if "Module enrichment" in rail else []) + ["Metabolite enzymes", "Pathway network", "Parameters used"]
    check(rail == expected, f"result panes (got {rail})")

    # ----------------------------------------------------------- overview
    print("=== overview ===")
    ov = page.inner_text(".multiomics-net-result")
    check(all(LABEL in ov for LABEL in ["Genomics", "Transcriptomics", "Proteomics", "Metabolomics"]) and "What each omics contributed" in ov, "the overview lists what every omics contributed")
    check(page.locator(".result-pane:visible tbody tr").count() == 4, "one row per omics in the contribution table")
    check(abs(sum(r["Share of the combined weight"] for r in R["omics"]) - 1.0) < 1e-9, "the shares of the combined weight sum to 1")
    with page.expect_download() as dl:
        page.locator("button", has_text="Download this table").click()
    rows = list(csv.DictReader(io.StringIO(open(dl.value.path(), encoding="utf-8-sig").read())))
    check(len(rows) == 4, "the contribution table downloads as a CSV with four omics")

    # ----------------------------------------------------------- hit table
    print("=== hit table ===")
    page.locator("button.result-nav-btn", has_text="Hit table").click()
    page.wait_for_timeout(700)
    n_rec = R["propagation"]["n_recruited"]
    cb = page.locator("input.multiomics-net-recruited-only")
    check(cb.is_checked() == (n_rec > 0), f"the table opens on the RECRUITED view exactly when there are recruited proteins ({n_rec})")
    if cb.is_checked():
        cb.click()
        page.wait_for_timeout(500)
    shown = page.locator(".result-pane:visible tbody tr").count()
    check(shown == min(len(R["propagation"]["table"]["rows"]), 2000), f"all proteins are listed (the first {shown} of {R['propagation']['n_proteins']}; the CSV has them all)")
    cols = [t.strip() for t in page.locator(".result-pane:visible thead th").all_inner_texts()]
    check("Seed omics" in cols and "Transcriptomics log2FC" in cols and "Proteomics padj" in cols, "the table names the seeding omics and each omics' own measurement")
    with page.expect_download() as dl:
        page.locator("button", has_text="Download hit table").click()
    n_csv = len(list(csv.DictReader(io.StringIO(open(dl.value.path(), encoding="utf-8-sig").read()))))
    check(n_csv == R["propagation"]["n_proteins"], f"the hit-table CSV has every protein ({n_csv})")

    # ----------------------------------------------------------- figures
    print("=== figures (the single-omics tiles' own builders) ===")
    open_pane(page, "Figures")
    n_links = len(R["propagation"]["edges"])
    page.locator("button", has_text="Show network figure").click()
    page.wait_for_timeout(2500)
    if n_links:
        check(page.locator(".result-pane:visible .js-plotly-plot").count() == 1, f"the propagation network figure draws from this payload ({n_links} STRING links among the significant proteins)")
    else:
        check(page.locator(".result-pane:visible .js-plotly-plot").count() == 0 and re.search(r"Not enough significant hits|No STRING links", page.inner_text(".result-pane:visible")),
              f"with {R['propagation']['n_hits']} significant proteins and no STRING link among them the figure SAYS so instead of drawing nothing")
    page.locator("button", has_text="Show null distribution").click()
    page.wait_for_timeout(2000)
    check(page.locator(".result-pane:visible .js-plotly-plot").count() == (2 if n_links else 1), "the null-distribution histogram draws")

    # ----------------------------------------------------------- modules
    mods = R["modules"]
    if mods and mods["table"]["rows"]:
        print("=== modules ===")
        open_pane(page, "Modules")
        mt_cols = [t.strip() for t in page.locator(".result-pane:visible thead th").all_inner_texts()]
        check(all(f"Seeds from {l}" in mt_cols for l in ["Genomics", "Transcriptomics", "Proteomics", "Metabolomics"]) and "Seeds from several omics" in mt_cols, "the module table counts the seeds of each omics")
        check(page.evaluate(HAS_DATA, "multiomics-net-bars"), "the module composition chart is drawn")
        names = page.evaluate("document.getElementById('multiomics-net-bars').data.map(t => t.name)")
        check(len(names) == 5 and names[-1].startswith("Recruited"), f"one bar series per omics plus the recruited share ({names})")
        check(page.evaluate(MODEBAR, "multiomics-net-bars"), "the composition chart has the figure tools row")
    else:
        print("   (no module at this run size - the module checks are skipped)")

    # ----------------------------------------------------------- metabolites
    print("=== metabolite enzymes ===")
    open_pane(page, "Metabolite enzymes")
    check(page.locator(".result-pane:visible tbody tr").count() == len(R["metabolites"]["rows"]), f"one row per hit metabolite ({len(R['metabolites']['rows'])})")
    txt = page.inner_text(".result-pane:visible")
    check("left out: in" in txt and "Family reactions skipped" in txt, "hubs are listed with the reason, and the family-reaction column is shown")

    # ----------------------------------------------------------- pathway network
    print("=== pathway network ===")
    open_pane(page, "Pathway network")
    pw = R["pathway_network"]
    check(page.evaluate(HAS_DATA, "multiomics-net-pathway"), "the pathway network is drawn")
    m = page.evaluate("""() => { const d = document.getElementById('multiomics-net-pathway'); const xa = d._fullLayout.xaxis, ya = d._fullLayout.yaxis;
        const nodeT = d.data.find(t => t.mode === 'markers+text');
        const pts = nodeT.x.map((x, i) => [xa.l2p(x), ya.l2p(nodeT.y[i])]);
        let minD = 1e9;
        for (let i = 0; i < pts.length; i++) for (let j = i + 1; j < pts.length; j++) minD = Math.min(minD, Math.hypot(pts[i][0] - pts[j][0], pts[i][1] - pts[j][1]));
        const labels = [...d.querySelectorAll('.textpoint text')].map(t => { const r = t.getBoundingClientRect(); return [r.left, r.top, r.right, r.bottom]; });
        let overlaps = 0;
        for (let i = 0; i < labels.length; i++) for (let j = i + 1; j < labels.length; j++) {
            const a = labels[i], b = labels[j];
            if (a[0] < b[2] - 1 && b[0] < a[2] - 1 && a[1] < b[3] - 1 && b[1] < a[3] - 1) overlaps++;
        }
        const chips = d.data.filter(t => t.mode === 'markers' && t.name).map(t => t.name);
        return {n: nodeT.x.length, minD, overlaps, nLabels: labels.length, chips, width: d.clientWidth, tool: !!d.querySelector('.modebar'),
                frame: d.parentElement.clientWidth, scroll: getComputedStyle(d.parentElement).overflow}; }""")
    check(m["n"] == len(pw["nodes"]), f"every pathway is drawn ({m['n']})")
    check(m["minD"] >= 40, f"no two pathways are closer than a marker plus its chips (closest {m['minD']:.0f}px)")
    check(m["nLabels"] == m["n"] and m["overlaps"] == 0, f"NO TWO LABELS OVERLAP (measured {m['nLabels']} label boxes on the page; {m['overlaps']} overlapping pairs)")
    check(sorted(m["chips"]) == sorted([{"genomics": "Genomics", "transcriptomics": "Transcriptomics", "proteomics": "Proteomics", "metabolomics": "Metabolomics"}[o] for o in pw["omics"]]),
          f"one evidence chip series per omics ({m['chips']})")
    check(m["tool"], "the pathway network has the figure tools row")
    check(m["scroll"] == "auto" and m["width"] >= 1100, f"drawn at a fixed readable size ({m['width']}px) inside a scroll frame ({m['frame']}px wide)")
    with page.expect_download() as dl:
        page.locator("button", has_text="Download pathways").click()
    rows = list(csv.DictReader(io.StringIO(open(dl.value.path(), encoding="utf-8-sig").read())))
    check(len(rows) == len(pw["nodes"]) and "Similar pathways merged into this one" in rows[0], f"the pathways CSV has every pathway and the merged ones ({len(rows)})")
    check(not pw["n_merged"] or f"{pw['n_merged']} further pathway(s) were merged" in page.inner_text(".result-pane:visible"), f"the pane says how many near-identical pathways were merged ({pw['n_merged']})")

    # ----------------------------------------------------------- navigation
    print("=== navigation ===")
    page.locator("button", has_text="Change settings").click()
    page.wait_for_timeout(400)
    check(page.evaluate("STATE.aaPage") == "multiomics_net" and page.evaluate("!!STATE.multiomics.resolved"), "Change settings returns to the setup with everything resolved kept")
    page.locator("button", has_text="Back to Multi-omics").click()
    page.wait_for_timeout(300)
    tiles = page.locator(".aa-genomics-grid button.module-btn")
    tiles.nth(3).click()
    page.wait_for_timeout(500)
    check(page.evaluate("STATE.aaPage") == "multiomics_net", "re-opening the box reuses the resolved uploads (no second upload)")

    print("=== console ===")
    check(len(errors) == 0, f"no page error or console error during the whole flow (got {errors})")
    browser.close()

print(f"\nTOTAL: {PASS} passed, {FAIL} failed")
if FAIL:
    raise SystemExit(1)
