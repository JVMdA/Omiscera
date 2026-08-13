"""Regression test for a real, user-reported bug in buildGseaNullHistogram
(static/app.js): when the observed NES is more extreme than the ENTIRE
permutation null distribution's own range - exactly what happens for a
strongly significant result - Plotly's default autorange stretches the
x-axis JUST far enough to include the observed-NES shape line, with
near-zero padding, pinning it directly onto the axis border where it is
visually indistinguishable from it. A real screenshot (NES=-2.40,
p=6.8e-3, a chromosome-organization GSEA term) showed no visible line at
all, even though the shapes array genuinely contained one.

THE FIX: buildGseaNullHistogram now computes an explicit xaxis.range that
pads out from whichever is more extreme - the null samples, the observed
NES itself, or 0 (the always-drawn zero-reference line) - by 6% of the
total span, guaranteeing real, visible clearance no matter how extreme the
observed NES is.

Needs a live server (to load app.js and Plotly) - no network calls,
synthetic fixtures.
"""
import random
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


def clearance_fraction(page, samples, nes):
    """Build the chart with the REAL buildGseaNullHistogram and return how
    far (as a fraction of the rendered axis span) the observed-NES shape
    line sits from the nearer axis edge - the direct measure of "is it
    visibly clear of the border" that the bug and its fix are both about."""
    page.evaluate(
        """([samples, nes]) => {
        window.__testDiv = document.createElement('div');
        window.__testDiv.style.width = '1100px';
        document.body.appendChild(window.__testDiv);
        const block = { null_nes_by_size: { '107': samples } };
        const termRow = { 'Term Name': 'test term', 'Size': 107, 'NES': nes, 'p-value': 0.0068 };
        buildGseaNullHistogram(window.__testDiv, block, termRow, 12, 'test');
    }""",
        [samples, nes],
    )
    page.wait_for_timeout(500)
    result = page.evaluate("""() => {
        const gd = window.__testDiv;
        return { xrange: gd._fullLayout.xaxis.range, shapeX: gd.layout.shapes[1].x0 };
    }""")
    page.evaluate("window.__testDiv.remove()")
    lo, hi = result["xrange"]
    span = hi - lo
    dist_from_edge = min(result["shapeX"] - lo, hi - result["shapeX"])
    return dist_from_edge / span if span else 0


with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page()
    errors = []
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    page.goto(BASE, wait_until="networkidle")
    page.wait_for_timeout(300)

    random.seed(0)
    samples = [random.gauss(-1, 0.25) for _ in range(250)] + [random.gauss(1, 0.2) for _ in range(250)]

    print("=== observed NES more extreme than the whole null distribution (the bug's exact trigger) ===")
    frac_extreme = clearance_fraction(page, samples, -2.40)
    check(frac_extreme > 0.02,
          f"the observed-NES line sits with real (>2%) clearance from the axis edge, not pinned to the border "
          f"(got {frac_extreme:.4f} of the axis span - the real pre-fix bug measured ~0.001)")

    print("=== observed NES well within the null distribution's own range (the common case) ===")
    frac_typical = clearance_fraction(page, samples, 0.5)
    check(frac_typical > 0.02,
          f"a typical, non-extreme observed NES also keeps real clearance (got {frac_typical:.4f})")

    print("=== SENSITIVITY: the clearance check is not vacuous ===")
    # Re-run the EXTREME case exactly as the PRE-FIX buildGseaNullHistogram
    # rendered it: same trace and shapes, but no explicit xaxis.range - left
    # entirely to Plotly's own autorange, which is what produced the
    # invisible line in the first place.
    page.evaluate(
        """([samples, nes]) => {
        window.__testDiv = document.createElement('div');
        window.__testDiv.style.width = '1100px';
        document.body.appendChild(window.__testDiv);
        const col = nes >= 0 ? '#B84A52' : '#1B639C';
        const trace = { x: samples, type: 'histogram', nbinsx: 60 };
        const layout = {
            height: 460,
            xaxis: { title: { text: 'x' }, zeroline: true },
            yaxis: { title: { text: 'y' } },
            shapes: [
                { type: 'line', x0: 0, x1: 0, yref: 'paper', y0: 0, y1: 1 },
                { type: 'line', x0: nes, x1: nes, yref: 'paper', y0: 0, y1: 1, line: { color: col, width: 2 } },
            ],
        };
        Plotly.newPlot(window.__testDiv, [trace], layout, {});
    }""",
        [samples, -2.40],
    )
    page.wait_for_timeout(500)
    broken_result = page.evaluate("""() => {
        const gd = window.__testDiv;
        return { xrange: gd._fullLayout.xaxis.range, shapeX: gd.layout.shapes[1].x0 };
    }""")
    page.evaluate("window.__testDiv.remove()")
    lo, hi = broken_result["xrange"]
    span = hi - lo
    dist_from_edge = min(broken_result["shapeX"] - lo, hi - broken_result["shapeX"])
    broken_frac = dist_from_edge / span if span else 0
    check(broken_frac < 0.02,
          f"the PRE-FIX rendering (no explicit range, Plotly's own autorange only) genuinely pins the line "
          f"within 2% of the axis edge for this exact input (got {broken_frac:.4f}) - proving the check above "
          f"is not vacuous")

    check(len(errors) == 0, f"no page errors while exercising buildGseaNullHistogram (got {errors})")
    browser.close()

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
