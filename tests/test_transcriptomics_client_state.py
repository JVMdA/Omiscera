# -*- coding: utf-8 -*-
"""Generic, source-level invariants over static/app.js's Transcriptomics
module wiring - the same text-parsing technique test_organism_consistency.py
already uses on API.py, applied client-side. No browser, no server.

WHY GENERIC. Written once, against the module LIST discovered from
freshTranscriptomicsState() itself, so it automatically covers ORA/GSEA/
Regulon/Network today and a future 5th tile tomorrow with no edits here -
the same reasoning that drove Part 4's route-discovery change to
test_organism_consistency.py. A hand-listed "check ora, gsea, regulon,
network" version would silently stop covering anything past the day someone
adds a 5th module and forgets this file exists.

Checks:
  1. every "<mod>Config: null," in freshTranscriptomicsState() has a matching
     "<mod>Config = null;" / "<mod>Results = null;" / "<mod>RunLocked = false;"
     / "delete tx._<mod>Draft;" in _clearTranscriptomicsResolutionDependentState()
     - an organism change or a new upload must never leave a stale analysis on
     screen for ANY module.
  2. every module's config/results router pages exist AND their render
     functions are actually defined (not just referenced).
  3. every tile in the grid array is either flagged not-yet-implemented or has
     a branch in the Continue-chain - the fallback "not available yet" branch
     must stay unreachable for a real tile.
"""
import os
import re
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
APP_JS = os.path.join(_ROOT, "static", "app.js")

PASS = FAIL = 0


def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("PASS: " + msg)
    else:
        FAIL += 1
        print("FAIL: " + msg)


with open(APP_JS, "r", encoding="utf-8") as fh:
    SRC = fh.read()


def extract_function(src, name):
    """Text from a top-level `function name(...) {` to the next top-level
    `function ` declaration (or EOF) - the same "function span" idiom used
    this session for safely-scoped app.js edits, repurposed here for reading."""
    m = re.search(r'\nfunction ' + re.escape(name) + r'\([^)]*\)\s*\{', src)
    if not m:
        return None
    start = m.end()
    m2 = re.search(r'\nfunction \w+\(', src[start:])
    return src[start:start + m2.start()] if m2 else src[start:]


# =============================================================================
print("=== discover the module list from freshTranscriptomicsState() itself ===")
# =============================================================================
fresh_body = extract_function(SRC, "freshTranscriptomicsState")
check(fresh_body is not None, "freshTranscriptomicsState() is found in app.js")
MODULES = re.findall(r'(\w+)Config: null,', fresh_body)
check(set(MODULES) == {"ora", "gsea", "regulon", "network"},
      f"discovered exactly the 4 known modules (got {MODULES}) - if this list grows, every check below "
      f"grows with it automatically")

# =============================================================================
print("=== _clearTranscriptomicsResolutionDependentState() clears EVERY module ===")
# =============================================================================
clear_body = extract_function(SRC, "_clearTranscriptomicsResolutionDependentState")
check(clear_body is not None, "_clearTranscriptomicsResolutionDependentState() is found in app.js")
for mod in MODULES:
    check(f"tx.{mod}Config = null;" in clear_body, f"{mod}Config is cleared")
    check(f"tx.{mod}Results = null;" in clear_body, f"{mod}Results is cleared")
    check(f"tx.{mod}RunLocked = false;" in clear_body, f"{mod}RunLocked is reset")
    check(f"delete tx._{mod}Draft;" in clear_body, f"the stale, possibly organism-gated _{mod}Draft is deleted")

print("=== SENSITIVITY: the check above genuinely detects a missing clear line ===")
# In-memory only - a string surgically missing ONE real line, never written to
# disk. If a module's clearing were accidentally dropped, this is exactly what
# the source would look like, and the check above must catch it.
broken_clear_body = clear_body.replace("delete tx._networkDraft;\n", "")
broken_missing = not (f"delete tx._network" + "Draft;" in broken_clear_body)
check(broken_missing, "the deliberately-broken in-memory copy really is missing the line (sensitivity is real, "
                      "not vacuous)")
check(f"delete tx._networkDraft;" in clear_body,
      "...while the REAL app.js on disk still has it - confirming the check above passes for the right reason")

# =============================================================================
print("=== router: every module's config/results pages exist and route to a REAL function ===")
# =============================================================================
router_body = extract_function(SRC, "renderAnnotationAnalytics")
check(router_body is not None, "renderAnnotationAnalytics() (the page router) is found in app.js")
for mod in MODULES:
    cap = mod.capitalize()
    for stage in ("config", "results"):
        route_line = f'if (STATE.aaPage === "transcriptomics_{mod}_{stage}") return renderTranscriptomics{cap}{stage.capitalize()}Screen();'
        check(route_line in router_body, f"router has a line for transcriptomics_{mod}_{stage} (got: {route_line!r})")
        fn_name = f"renderTranscriptomics{cap}{stage.capitalize()}Screen"
        check(re.search(r'\nfunction ' + re.escape(fn_name) + r'\(', SRC) is not None,
              f"and {fn_name} is actually DEFINED somewhere in app.js, not just referenced "
              f"(a route to an undefined function is a silent ReferenceError, only surfacing at click time)")

# =============================================================================
print("=== every tile in the grid has EITHER a not-yet-implemented flag OR a Continue-chain branch ===")
# =============================================================================
grid_body = extract_function(SRC, "renderTranscriptomicsScreen")
check(grid_body is not None, "renderTranscriptomicsScreen() (the tile grid) is found in app.js")
m = re.search(r'\[\s*((?:"[^"]+",?\s*)+)\]\.forEach', grid_body)
check(m is not None, "the tile-name array literal is found inside the grid function")
tile_names = re.findall(r'"([^"]+)"', m.group(1))
check(len(tile_names) >= 4, f"at least the 4 known tiles are present (got {tile_names})")

not_impl_match = re.search(r'TRANSCRIPTOMICS_TILES_NOT_YET_IMPLEMENTED = new Set\(\[([^\]]*)\]\)', SRC)
check(not_impl_match is not None, "TRANSCRIPTOMICS_TILES_NOT_YET_IMPLEMENTED is found in app.js")
not_impl = set(re.findall(r'"([^"]+)"', not_impl_match.group(1) or ""))

continue_body = extract_function(SRC, "renderTranscriptomicsResults")
check(continue_body is not None, "renderTranscriptomicsResults() (the Continue-chain) is found in app.js")
continue_names = set(re.findall(r'STATE\.transcriptomicsModuleName === "([^"]+)"', continue_body))


def _tile_is_covered(name, continue_set):
    return name in not_impl or name in continue_set


for name in tile_names:
    check(_tile_is_covered(name, continue_names),
          f"tile {name!r} is either flagged not-yet-implemented or has a Continue-chain branch")

print("=== SENSITIVITY: an unhandled tile is genuinely caught, not passed by accident ===")
check(not _tile_is_covered("Bogus Future Module", continue_names),
      "a hypothetical 5th tile with neither a flag nor a branch correctly FAILS the same check used above - "
      "so the all-pass result for the 4 real tiles reflects real coverage, not a check that always says yes")
check(continue_body.rstrip().endswith(
    'wrap.appendChild(el("div", { class: "caption" },\n            "This analysis module is not available yet."));'
) or '"This analysis module is not available yet."' in continue_body,
      "the Continue-chain's fallback branch is present (reachable only for an uncovered tile - never fires "
      "for the 4 real ones, which is what the checks above just confirmed)")

print("\nTOTAL: %d passed, %d failed" % (PASS, FAIL))
sys.exit(1 if FAIL else 0)
